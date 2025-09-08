#ifndef HDF5_UTILS_HPP
#define HDF5_UTILS_HPP

#include <Eigen/Dense>
#include <string>
#include <H5Cpp.h>
#include <iostream>

/**
 * HDF5 utility class for saving Eigen matrices and vectors
 * Provides a clean interface for saving simulation results to HDF5 format
 */
class HDF5Utils {
public:
    /**
     * Save a single Eigen matrix to HDF5 file
     * 
     * @param matrix: Eigen matrix to save
     * @param filename: output HDF5 file name
     * @param dataset_name: name of the dataset in HDF5 file
     * @param overwrite: whether to overwrite existing file (default: true)
     */
    static void saveMatrix(
        const Eigen::MatrixXd& matrix,
        const std::string& filename,
        const std::string& dataset_name,
        bool overwrite = true
    );

    /**
     * Save a single Eigen vector to HDF5 file
     * 
     * @param vector: Eigen vector to save
     * @param filename: output HDF5 file name
     * @param dataset_name: name of the dataset in HDF5 file
     * @param overwrite: whether to overwrite existing file (default: true)
     */
    static void saveVector(
        const Eigen::VectorXd& vector,
        const std::string& filename,
        const std::string& dataset_name,
        bool overwrite = true
    );

    /**
     * Save travel time results with coordinates
     * 
     * @param xx: x coordinates vector
     * @param yy: y coordinates vector
     * @param T: travel time matrix
     * @param filename: output HDF5 file name
     */
    static void saveTravelTimeResults(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& yy,
        const Eigen::MatrixXd& T,
        const std::string& filename
    );

    /**
     * Save complete simulation results including coefficients
     * 
     * @param xx: x coordinates vector
     * @param yy: y coordinates vector
     * @param T: travel time matrix
     * @param filename: output HDF5 file name
     * @param a: coefficient matrix a (optional)
     * @param b: coefficient matrix b (optional)
     * @param c: coefficient matrix c (optional)
     * @param slowness: slowness field (optional)
     * @param source_x: source x coordinate (optional)
     * @param source_y: source y coordinate (optional)
     */
    static void saveCompleteResults(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& yy,
        const Eigen::MatrixXd& T,
        const std::string& filename,
        const Eigen::MatrixXd* a = nullptr,
        const Eigen::MatrixXd* b = nullptr,
        const Eigen::MatrixXd* c = nullptr,
        const Eigen::MatrixXd* slowness = nullptr,
        const double* source_x = nullptr,
        const double* source_y = nullptr
    );

    /**
     * Add attributes to HDF5 file (metadata)
     * 
     * @param filename: HDF5 file name
     * @param attr_name: attribute name
     * @param attr_value: attribute value (string)
     */
    static void addStringAttribute(
        const std::string& filename,
        const std::string& attr_name,
        const std::string& attr_value
    );

    /**
     * Add numeric attribute to HDF5 file
     * 
     * @param filename: HDF5 file name
     * @param attr_name: attribute name
     * @param attr_value: attribute value (double)
     */
    static void addDoubleAttribute(
        const std::string& filename,
        const std::string& attr_name,
        double attr_value
    );

private:
    /**
     * Helper function to create or open HDF5 file
     */
    static H5::H5File createOrOpenFile(const std::string& filename, bool overwrite = true);

    /**
     * Helper function to save matrix data to existing file
     */
    static void saveMatrixToFile(
        H5::H5File& file,
        const Eigen::MatrixXd& matrix,
        const std::string& dataset_name
    );

    /**
     * Helper function to save vector data to existing file
     */
    static void saveVectorToFile(
        H5::H5File& file,
        const Eigen::VectorXd& vector,
        const std::string& dataset_name
    );
};

#endif // HDF5_UTILS_HPP
