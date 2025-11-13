//****h* lut_mpi_hdf5.cc
// NAME
//   lut_mpi_hdf5.cc - MPI-based LUT generation with HDF5 output
//
// DESCRIPTION
//   MPI parallel LUT generation using MODTRAN C++ API with direct HDF5 output.
//   Master process distributes parameter combinations to workers, gathers results,
//   and writes consolidated HDF5 file.
//
// USAGE
//   mpirun -np <nproc> ./lut_mpi_hdf5 template.json output.h5 albedo_min albedo_max n_albedo h2o_min h2o_max n_h2o
//
// EXAMPLE
//   mpirun -np 4 ./lut_mpi_hdf5 modtran_template.json lut.h5 0.0 1.0 11 0.0 4.0 151
//****

#include <modtran/modtran.h>
#include <nlohmann/json.hpp>
#include <hdf5.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>

using json = nlohmann::json;

#define MASTER (0)
#define TAG_TASK (1)
#define TAG_RESULT (2)
#define TAG_DONE (999)

// Task structure sent from master to workers
struct Task {
  int task_id;
  int albedo_idx;
  int h2o_idx;
  double albedo;
  double h2o;
};

// Result structure sent from workers to master
struct ResultHeader {
  int task_id;
  int albedo_idx;
  int h2o_idx;
  int num_channels;
  int success;
};

// Read JSON file
std::string readJsonFile(const char* filepath)
{
  std::ifstream file(filepath);
  if (!file.is_open())
  {
    std::cerr << "Error: Could not open file " << filepath << std::endl;
    return "";
  }
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

// Generate linear space
std::vector<double> linspace(double start, double end, int num)
{
  std::vector<double> result(num);
  if (num == 1) {
    result[0] = start;
  } else {
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++)
      result[i] = start + i * step;
  }
  return result;
}

// HDF5 helper: Write 1D array
bool writeHDF5Array(hid_t file_id, const char* name, const std::vector<double>& data)
{
  hsize_t dims[1] = {(hsize_t)data.size()};
  hid_t dataspace = H5Screate_simple(1, dims, NULL);
  hid_t dataset = H5Dcreate2(file_id, name, H5T_NATIVE_DOUBLE, dataspace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
    H5Sclose(dataspace);
    return false;
  }
  herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, data.data());
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return (status >= 0);
}

// HDF5 helper: Write 3D array [albedo][h2o][channel]
bool writeHDF5Dataset(hid_t file_id, const char* name,
                      const std::vector<std::vector<std::vector<float>>>& data,
                      int n_albedo, int n_h2o, int n_chan)
{
  // Flatten 3D to 1D
  std::vector<float> flat_data;
  flat_data.reserve(n_albedo * n_h2o * n_chan);

  for (int i = 0; i < n_albedo; i++)
    for (int j = 0; j < n_h2o; j++)
      for (int k = 0; k < n_chan; k++)
        flat_data.push_back(data[i][j][k]);

  hsize_t dims[3] = {(hsize_t)n_albedo, (hsize_t)n_h2o, (hsize_t)n_chan};
  hid_t dataspace = H5Screate_simple(3, dims, NULL);
  hid_t dataset = H5Dcreate2(file_id, name, H5T_NATIVE_FLOAT, dataspace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
    H5Sclose(dataspace);
    return false;
  }
  herr_t status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, flat_data.data());
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return (status >= 0);
}

int main(int argc, char* argv[])
{
  int my_rank, nproc;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (argc < 9)
  {
    if (my_rank == MASTER)
    {
      std::cout << "Usage: " << argv[0] << " template.json output.h5 "
                << "albedo_min albedo_max n_albedo h2o_min h2o_max n_h2o" << std::endl;
      std::cout << "Example: mpirun -np 4 " << argv[0] << " template.json lut.h5 0.0 1.0 11 0.0 4.0 151" << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  const char* template_path = argv[1];
  const char* output_path = argv[2];
  double albedo_min = atof(argv[3]);
  double albedo_max = atof(argv[4]);
  int n_albedo = atoi(argv[5]);
  double h2o_min = atof(argv[6]);
  double h2o_max = atof(argv[7]);
  int n_h2o = atoi(argv[8]);

  if (nproc < 2)
  {
    if (my_rank == MASTER)
      std::cerr << "ERROR: At least 2 processes required (1 master + 1 worker)" << std::endl;
    MPI_Finalize();
    return -1;
  }

  // Load JSON template (all processes)
  std::string json_str = readJsonFile(template_path);
  if (json_str.empty())
  {
    MPI_Finalize();
    return -1;
  }
  json template_config = json::parse(json_str);

  // Generate parameter grids
  std::vector<double> albedo_grid = linspace(albedo_min, albedo_max, n_albedo);
  std::vector<double> h2o_grid = linspace(h2o_min, h2o_max, n_h2o);
  int total_tasks = n_albedo * n_h2o;

  if (my_rank == MASTER)
  {
    std::cout << "\n========================================" << std::endl;
    std::cout << "MPI + HDF5 LUT Generation" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nTemplate:  " << template_path << std::endl;
    std::cout << "Output:    " << output_path << std::endl;
    std::cout << "Processes: " << nproc << " (" << (nproc-1) << " workers)" << std::endl;
    std::cout << "\nParameter grid: " << n_albedo << " albedos × " << n_h2o << " H2O = "
              << total_tasks << " runs" << std::endl;
    std::cout << "Albedo: [" << albedo_min << ", " << albedo_max << "]" << std::endl;
    std::cout << "H2O:    [" << h2o_min << ", " << h2o_max << "] g/cm²" << std::endl;

    // Allocate result storage (assume 425 channels for AVIRIS-NG)
    int num_channels = 425;
    std::vector<std::vector<std::vector<float>>> total_radiance(n_albedo,
      std::vector<std::vector<float>>(n_h2o, std::vector<float>(num_channels, 0.0f)));
    std::vector<std::vector<std::vector<float>>> transmission(n_albedo,
      std::vector<std::vector<float>>(n_h2o, std::vector<float>(num_channels, 0.0f)));
    std::vector<std::vector<std::vector<float>>> channel_rad(n_albedo,
      std::vector<std::vector<float>>(n_h2o, std::vector<float>(num_channels, 0.0f)));
    std::vector<std::vector<std::vector<float>>> grnd_rflt(n_albedo,
      std::vector<std::vector<float>>(n_h2o, std::vector<float>(num_channels, 0.0f)));
    std::vector<std::vector<std::vector<float>>> mult_scat(n_albedo,
      std::vector<std::vector<float>>(n_h2o, std::vector<float>(num_channels, 0.0f)));
    std::vector<std::vector<std::vector<float>>> sing_scat(n_albedo,
      std::vector<std::vector<float>>(n_h2o, std::vector<float>(num_channels, 0.0f)));

    std::cout << "\nAllocated storage for " << n_albedo << " × " << n_h2o << " × "
              << num_channels << " channels" << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "DISTRIBUTING TASKS TO WORKERS" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Distribute tasks
    int next_task = 0;
    int completed = 0;
    MPI_Status stat;

    // Initially send one task to each worker
    for (int worker = 1; worker < nproc && next_task < total_tasks; worker++)
    {
      int i = next_task / n_h2o;
      int j = next_task % n_h2o;
      Task task = {next_task, i, j, albedo_grid[i], h2o_grid[j]};
      MPI_Send(&task, sizeof(Task), MPI_BYTE, worker, TAG_TASK, MPI_COMM_WORLD);
      std::cout << "Task " << next_task << " -> Worker " << worker
                << " (albedo=" << std::fixed << std::setprecision(3) << task.albedo
                << ", H2O=" << task.h2o << ")" << std::endl;
      next_task++;
    }

    // Receive results and send more tasks
    while (completed < total_tasks)
    {
      // Receive result header
      ResultHeader header;
      MPI_Recv(&header, sizeof(ResultHeader), MPI_BYTE, MPI_ANY_SOURCE, TAG_RESULT,
               MPI_COMM_WORLD, &stat);
      int worker = stat.MPI_SOURCE;

      if (header.success)
      {
        // Receive channel data
        int i = header.albedo_idx;
        int j = header.h2o_idx;

        MPI_Recv(total_radiance[i][j].data(), num_channels, MPI_FLOAT, worker, TAG_RESULT, MPI_COMM_WORLD, &stat);
        MPI_Recv(transmission[i][j].data(), num_channels, MPI_FLOAT, worker, TAG_RESULT, MPI_COMM_WORLD, &stat);
        MPI_Recv(channel_rad[i][j].data(), num_channels, MPI_FLOAT, worker, TAG_RESULT, MPI_COMM_WORLD, &stat);
        MPI_Recv(grnd_rflt[i][j].data(), num_channels, MPI_FLOAT, worker, TAG_RESULT, MPI_COMM_WORLD, &stat);
        MPI_Recv(mult_scat[i][j].data(), num_channels, MPI_FLOAT, worker, TAG_RESULT, MPI_COMM_WORLD, &stat);
        MPI_Recv(sing_scat[i][j].data(), num_channels, MPI_FLOAT, worker, TAG_RESULT, MPI_COMM_WORLD, &stat);

        completed++;
        std::cout << "Task " << header.task_id << " completed by Worker " << worker
                  << " (" << completed << "/" << total_tasks << " = "
                  << std::fixed << std::setprecision(1) << (100.0 * completed / total_tasks) << "%)" << std::endl;
      }
      else
      {
        std::cerr << "Task " << header.task_id << " FAILED on Worker " << worker << std::endl;
        completed++;  // Count as done to avoid hanging
      }

      // Send next task if available
      if (next_task < total_tasks)
      {
        int i = next_task / n_h2o;
        int j = next_task % n_h2o;
        Task task = {next_task, i, j, albedo_grid[i], h2o_grid[j]};
        MPI_Send(&task, sizeof(Task), MPI_BYTE, worker, TAG_TASK, MPI_COMM_WORLD);
        std::cout << "Task " << next_task << " -> Worker " << worker
                  << " (albedo=" << std::fixed << std::setprecision(3) << task.albedo
                  << ", H2O=" << task.h2o << ")" << std::endl;
        next_task++;
      }
      else
      {
        // No more tasks, send done signal
        Task done_task = {-1, 0, 0, 0.0, 0.0};
        MPI_Send(&done_task, sizeof(Task), MPI_BYTE, worker, TAG_DONE, MPI_COMM_WORLD);
      }
    }

    // Send done signal to any remaining workers
    for (int worker = 1; worker < nproc; worker++)
    {
      Task done_task = {-1, 0, 0, 0.0, 0.0};
      MPI_Send(&done_task, sizeof(Task), MPI_BYTE, worker, TAG_DONE, MPI_COMM_WORLD);
    }

    // Write HDF5 file
    std::cout << "\n========================================" << std::endl;
    std::cout << "WRITING HDF5 OUTPUT" << std::endl;
    std::cout << "========================================\n" << std::endl;

    hid_t file_id = H5Fcreate(output_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
    {
      std::cerr << "ERROR: Could not create HDF5 file " << output_path << std::endl;
      MPI_Finalize();
      return -1;
    }

    writeHDF5Array(file_id, "albedo_grid", albedo_grid);
    writeHDF5Array(file_id, "h2o_grid", h2o_grid);
    std::cout << "  - Parameter grids" << std::endl;

    writeHDF5Dataset(file_id, "total_radiance", total_radiance, n_albedo, n_h2o, num_channels);
    std::cout << "  - total_radiance (" << n_albedo << " × " << n_h2o << " × " << num_channels << ")" << std::endl;

    writeHDF5Dataset(file_id, "transmission", transmission, n_albedo, n_h2o, num_channels);
    std::cout << "  - transmission" << std::endl;

    writeHDF5Dataset(file_id, "channel_radiance", channel_rad, n_albedo, n_h2o, num_channels);
    std::cout << "  - channel_radiance" << std::endl;

    writeHDF5Dataset(file_id, "grnd_rflt", grnd_rflt, n_albedo, n_h2o, num_channels);
    std::cout << "  - grnd_rflt" << std::endl;

    writeHDF5Dataset(file_id, "mult_scat", mult_scat, n_albedo, n_h2o, num_channels);
    std::cout << "  - mult_scat" << std::endl;

    writeHDF5Dataset(file_id, "sing_scat", sing_scat, n_albedo, n_h2o, num_channels);
    std::cout << "  - sing_scat" << std::endl;

    H5Fclose(file_id);

    std::cout << "\n========================================" << std::endl;
    std::cout << "SUCCESS!" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nCompleted: " << completed << " / " << total_tasks << " runs" << std::endl;
    std::cout << "HDF5 output: " << output_path << std::endl;
    std::cout << "\nVerify with Python:" << std::endl;
    std::cout << "  import h5py" << std::endl;
    std::cout << "  f = h5py.File('" << output_path << "', 'r')" << std::endl;
    std::cout << "  print(f.keys())" << std::endl;
    std::cout << "  print(f['total_radiance'].shape)  # Should be ("
              << n_albedo << ", " << n_h2o << ", " << num_channels << ")" << std::endl;
  }
  else  // WORKER processes
  {
    while (true)
    {
      Task task;
      MPI_Status stat;
      MPI_Recv(&task, sizeof(Task), MPI_BYTE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

      if (stat.MPI_TAG == TAG_DONE || task.task_id < 0)
        break;

      // Create thread-local JSON config
      json worker_config = template_config;

      // Generate unique NAME
      std::ostringstream name_stream;
      name_stream << "lut_" << std::setfill('0') << std::setw(2) << task.albedo_idx
                  << "_" << std::setw(2) << task.h2o_idx;
      worker_config["MODTRAN"][0]["MODTRANINPUT"]["NAME"] = name_stream.str();

      // Update parameters
      worker_config["MODTRAN"][0]["SURFACE"]["SURREF"] = task.albedo;
      worker_config["MODTRAN"][0]["ATMOSPHERE"]["H2OSTR"] = task.h2o;

      // Disable file outputs
      worker_config["MODTRAN"][0]["MODTRANINPUT"]["FILEOPTIONS"]["NOFILE"] = "FC_NOFILES";

      std::string config_str = worker_config.dump();

      // Execute MODTRAN
      modtran::ModLib modlib;
      int case_indx = modlib.inputJsonDoc(config_str.c_str(), config_str.length());

      ResultHeader header;
      header.task_id = task.task_id;
      header.albedo_idx = task.albedo_idx;
      header.h2o_idx = task.h2o_idx;
      header.num_channels = 425;
      header.success = 0;

      if (case_indx >= 0 && modlib.execute())
      {
        ModOutput* modout = modlib.caseOutput(case_indx);
        if (modout)
        {
          ModOutChanRadiance* chan = &(modout->channel_spectra.radiance);
          header.success = 1;

          // Send header
          MPI_Send(&header, sizeof(ResultHeader), MPI_BYTE, MASTER, TAG_RESULT, MPI_COMM_WORLD);

          // Send channel data
          MPI_Send(chan->total_rad_wavln, header.num_channels, MPI_FLOAT, MASTER, TAG_RESULT, MPI_COMM_WORLD);
          MPI_Send(chan->tot_trans, header.num_channels, MPI_FLOAT, MASTER, TAG_RESULT, MPI_COMM_WORLD);
          MPI_Send(chan->channel_rad, header.num_channels, MPI_FLOAT, MASTER, TAG_RESULT, MPI_COMM_WORLD);
          MPI_Send(chan->grnd_rflt, header.num_channels, MPI_FLOAT, MASTER, TAG_RESULT, MPI_COMM_WORLD);
          MPI_Send(chan->mult_scat, header.num_channels, MPI_FLOAT, MASTER, TAG_RESULT, MPI_COMM_WORLD);
          MPI_Send(chan->sing_scat, header.num_channels, MPI_FLOAT, MASTER, TAG_RESULT, MPI_COMM_WORLD);
        }
        else
        {
          MPI_Send(&header, sizeof(ResultHeader), MPI_BYTE, MASTER, TAG_RESULT, MPI_COMM_WORLD);
        }
      }
      else
      {
        MPI_Send(&header, sizeof(ResultHeader), MPI_BYTE, MASTER, TAG_RESULT, MPI_COMM_WORLD);
      }
    }
  }

  MPI_Finalize();
  return 0;
}
