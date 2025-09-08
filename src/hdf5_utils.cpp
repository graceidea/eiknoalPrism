#include "hdf5_utils.hpp"

void HDF5Utils::saveMatrix(
    const Eigen::MatrixXd& matrix,
    const std::string& filename,
    const std::string& dataset_name,
    bool overwrite
) {
    try {
        H5::H5File file = createOrOpenFile(filename, overwrite);
        saveMatrixToFile(file, matrix, dataset_name);
        std::cout << "Matrix saved to " << filename << " (dataset: " << dataset_name << ")" << std::endl;
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error while saving matrix to " << filename << ": " << e.getDetailMsg() << std::endl;
    }
}

void HDF5Utils::saveVector(
    const Eigen::VectorXd& vector,
    const std::string& filename,
    const std::string& dataset_name,
    bool overwrite
) {
    try {
        H5::H5File file = createOrOpenFile(filename, overwrite);
        saveVectorToFile(file, vector, dataset_name);
        std::cout << "Vector saved to " << filename << " (dataset: " << dataset_name << ")" << std::endl;
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error while saving vector to " << filename << ": " << e.getDetailMsg() << std::endl;
    }
}

void HDF5Utils::saveTravelTimeResults(
    const Eigen::VectorXd& xx,
    const Eigen::VectorXd& yy,
    const Eigen::MatrixXd& T,
    const std::string& filename
) {
    try {
        H5::H5File file = createOrOpenFile(filename, true);
        
        // Save coordinate vectors
        saveVectorToFile(file, xx, "xx");
        saveVectorToFile(file, yy, "yy");
        
        // Save travel time matrix
        saveMatrixToFile(file, T, "T");
        
        // Add metadata
        addStringAttribute(filename, "description", "Travel time simulation results");
        addStringAttribute(filename, "coordinates", "xx, yy are coordinate vectors");
        addStringAttribute(filename, "data", "T is the travel time matrix");
        
        std::cout << "Travel time results saved to " << filename << std::endl;
        
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error while saving travel time results to " << filename << ": " << e.getDetailMsg() << std::endl;
    }
}

void HDF5Utils::saveCompleteResults(
    const Eigen::VectorXd& xx,
    const Eigen::VectorXd& yy,
    const Eigen::MatrixXd& T,
    const std::string& filename,
    const Eigen::MatrixXd* a,
    const Eigen::MatrixXd* b,
    const Eigen::MatrixXd* c,
    const Eigen::MatrixXd* slowness,
    const double* source_x,
    const double* source_y
) {
    try {
        H5::H5File file = createOrOpenFile(filename, true);
        
        // Save coordinate vectors
        saveVectorToFile(file, xx, "xx");
        saveVectorToFile(file, yy, "yy");
        
        // Save travel time matrix
        saveMatrixToFile(file, T, "T");
        
        // Save optional coefficient matrices
        if (a != nullptr) {
            saveMatrixToFile(file, *a, "a");
        }
        if (b != nullptr) {
            saveMatrixToFile(file, *b, "b");
        }
        if (c != nullptr) {
            saveMatrixToFile(file, *c, "c");
        }
        if (slowness != nullptr) {
            saveMatrixToFile(file, *slowness, "slowness");
        }
        
        // Save source coordinates as attributes if provided
        if (source_x != nullptr) {
            addDoubleAttribute(filename, "source_x", *source_x);
        }
        if (source_y != nullptr) {
            addDoubleAttribute(filename, "source_y", *source_y);
        }
        
        // Add metadata
        addStringAttribute(filename, "description", "Complete Eikonal solver simulation results");
        addStringAttribute(filename, "solver", "FSM_WENO3_PS_2d");
        
        std::cout << "Complete results saved to " << filename << std::endl;
        
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error while saving complete results to " << filename << ": " << e.getDetailMsg() << std::endl;
    }
}

void HDF5Utils::addStringAttribute(
    const std::string& filename,
    const std::string& attr_name,
    const std::string& attr_value
) {
    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        
        // Create string datatype
        H5::StrType str_type(H5::PredType::C_S1, attr_value.size());
        
        // Create dataspace for attribute
        H5::DataSpace attr_space(H5S_SCALAR);
        
        // Create attribute
        H5::Attribute attr = file.createAttribute(attr_name, str_type, attr_space);
        
        // Write attribute
        attr.write(str_type, attr_value);
        
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error while adding string attribute " << attr_name << ": " << e.getDetailMsg() << std::endl;
    }
}

void HDF5Utils::addDoubleAttribute(
    const std::string& filename,
    const std::string& attr_name,
    double attr_value
) {
    try {
        H5::H5File file(filename, H5F_ACC_RDWR);
        
        // Create dataspace for attribute
        H5::DataSpace attr_space(H5S_SCALAR);
        
        // Create attribute
        H5::Attribute attr = file.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, attr_space);
        
        // Write attribute
        attr.write(H5::PredType::NATIVE_DOUBLE, &attr_value);
        
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 error while adding double attribute " << attr_name << ": " << e.getDetailMsg() << std::endl;
    }
}

H5::H5File HDF5Utils::createOrOpenFile(const std::string& filename, bool overwrite) {
    if (overwrite) {
        return H5::H5File(filename, H5F_ACC_TRUNC);
    } else {
        return H5::H5File(filename, H5F_ACC_RDWR);
    }
}

void HDF5Utils::saveMatrixToFile(
    H5::H5File& file,
    const Eigen::MatrixXd& matrix,
    const std::string& dataset_name
) {
    // Get matrix dimensions
    int nx = static_cast<int>(matrix.rows());
    int ny = static_cast<int>(matrix.cols());
    
    // Create dataspace with matrix dimensions
    hsize_t dims[2] = {static_cast<hsize_t>(nx), static_cast<hsize_t>(ny)};
    H5::DataSpace dataspace(2, dims);
    
    // Create dataset
    H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
    
    // Write data (Eigen matrices are column-major, HDF5 expects row-major)
    // We transpose the matrix to match HDF5 convention
    Eigen::MatrixXd matrix_transposed = matrix.transpose();
    dataset.write(matrix_transposed.data(), H5::PredType::NATIVE_DOUBLE);
}

void HDF5Utils::saveVectorToFile(
    H5::H5File& file,
    const Eigen::VectorXd& vector,
    const std::string& dataset_name
) {
    // Get vector dimensions
    hsize_t dims[1] = {static_cast<hsize_t>(vector.size())};
    H5::DataSpace dataspace(1, dims);
    
    // Create dataset
    H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
    
    // Write data
    dataset.write(vector.data(), H5::PredType::NATIVE_DOUBLE);
}
