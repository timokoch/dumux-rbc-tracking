// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup Microcirculation
 * \brief A red blood cell vtk writer
 * \author Timo Koch
 */

#ifndef DUMUX_MICROCIRCULATION_RBC_NETWORK_VTK_WRITER_HH
#define DUMUX_MICROCIRCULATION_RBC_NETWORK_VTK_WRITER_HH

#include <memory>
#include <fstream>
#include <iomanip>
#include <vector>

namespace Dumux {

/**
 * \file
 * \ingroup Microcirculation
 * \brief A vtk writer that writes out RBCs as cylinder segments for visualization
 * \author Timo Koch
 */
template <class RBCNetwork>
class RBCNetworkVTKWriter
{
    static constexpr unsigned int numBeforeLineBreak = 5;
public:

    RBCNetworkVTKWriter(std::shared_ptr<const RBCNetwork> rbcNetwork)
    : rbcNetwork_(rbcNetwork)
    {}

    RBCNetworkVTKWriter(std::shared_ptr<RBCNetwork> rbcNetwork)
    : rbcNetwork_(rbcNetwork)
    {}

    void write(const std::string& name) const
    {
        std::ofstream vtpFile(name + ".vtp");
        writeHeader_(vtpFile, rbcNetwork_->size(), rbcNetwork_->size()*2);

        // points
        writeCoordinates_(vtpFile);

        // lines
        writeLines_(vtpFile);

        // cell data
        vtpFile << "<CellData Scalars=\"radius\">\n";
        writeRadius_(vtpFile);
        for (int i = 0; i < cellData_.size(); ++i)
            writeCellData_(vtpFile, names_[i], *(cellData_[i]));
        vtpFile << "</CellData>\n";

        writeFooter_(vtpFile);
    }

    //! add cell data vector
    void addCellData(const std::vector<double>& cellData, const std::string& name)
    {
        names_.push_back(name);
        cellData_.push_back(&cellData);
    }

private:
    std::shared_ptr<const RBCNetwork> rbcNetwork_;

    std::vector<std::string> names_;
    std::vector<const std::vector<double>*> cellData_;

    /*!
     * \brief Writes the header to the file
     */
    void writeHeader_(std::ostream& file, std::size_t numLines, std::size_t numVertices) const
    {
        std::string header = "<?xml version=\"1.0\"?>\n";
                    header += "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                    header += "<PolyData>\n";
                    header += "<Piece NumberOfLines=\"" + std::to_string(numLines) + "\" NumberOfPoints=\"" + std::to_string(numVertices) + "\">\n";
        file << header;
    }

    /*!
     * \brief Writes the coordinates to the file
     * \param file The output file
     */
    void writeCoordinates_(std::ostream& file) const
    {
        // write the positions to the file
        file << "<Points>\n";
        file << "<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        int counter = 0;
        for (const auto& rbc : particles(*rbcNetwork_))
        {
            const auto length = rbcNetwork_->rbcLength(rbc);
            const auto& orientation = rbcNetwork_->orientation(rbc);
            auto pos = rbc.position();
            pos.axpy(-length*0.5, orientation);
            file << pos << " ";
            pos.axpy(length, orientation);
            file << pos << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
        file << "</Points>\n";
    }

    /*!
     * \brief Writes the lines to the file
     * \param file The output file
     */
    void writeLines_(std::ostream& file) const
    {
        // write the positions to the file
        file << "<Lines>\n";
        file << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int counter = 0;
        int vertexIdx = 0;
        for (std::size_t i = 0; i < rbcNetwork_->size(); ++i)
        {
            file << vertexIdx++ << " ";
            file << vertexIdx++ << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
        file << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        counter = 0;
        int vertexOffset = 0;
        for (std::size_t i = 0; i < rbcNetwork_->size(); ++i)
        {
            vertexOffset += 2;
            file << vertexOffset << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
        file << "</Lines>\n";
    }

    void writeRadius_(std::ostream& file) const
    {
        file << "<DataArray type=\"Float32\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int counter = 0;
        for (const auto& rbc : particles(*rbcNetwork_))
        {
            file << rbcNetwork_->radius(rbc) << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
    }

    void writeID_(std::ostream& file) const
    {
        file << "<DataArray type=\"Int64\" Name=\"id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int counter = 0;
        for (const auto& rbc : particles(*rbcNetwork_))
        {
            file << rbc.id() << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
    }

    void writeCellData_(std::ostream& file, const std::string& name, const std::vector<double>& data) const
    {
        file << "<DataArray type=\"Float32\" Name=\"" << name << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int counter = 0;
        for (const auto& rbc : particles(*rbcNetwork_))
        {
            file << data[rbc.id()] << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
    }

    void writeElementIdx_(std::ostream& file) const
    {
        file << "<DataArray type=\"Int64\" Name=\"eIdx\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int counter = 0;
        for (const auto& rbc : particles(*rbcNetwork_))
        {
            file << rbcNetwork_->elementIdx(rbc) << " ";

            // introduce a line break after a certain time
            if((++counter)  > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
    }

    /*!
     * \brief Writes the footer to the file
     */
    void writeFooter_(std::ostream& file) const
    {
        file << "</Piece>\n";
        file << "</PolyData>\n";
        file << "</VTKFile>";
    }
};

/**
 * \file
 * \ingroup Microcirculation
 * \brief A pvd writer (vtk time series) that writes out RBCs as cylinder segments for visualization
 * \author Timo Koch
 */
template <class RBCNetwork>
class RBCNetworkPVDWriter
{
public:
    RBCNetworkPVDWriter(std::shared_ptr<RBCNetwork> rbcNetwork, const std::string& name)
    : vtkWriter_(rbcNetwork)
    , name_(name)
    {}

    //! write a time point
    void write(double time)
    {
        const auto count = timeSteps_.size();
        timeSteps_.push_back(time);

        // write vtp file
        vtkWriter_.write(seqName_(count));

        // write pvd file
        std::ofstream pvdFile;
        pvdFile.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                           std::ios_base::eofbit);
        std::string pvdname = name_ + ".pvd";
        pvdFile.open(pvdname.c_str());
        pvdFile << "<?xml version=\"1.0\"?> \n"
                << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"> \n"
                << "<Collection> \n";

        for (std::size_t i = 0; i <= count; i++)
        {
            pvdFile << "<DataSet timestep=\"" << timeSteps_[i]
                    << "\" group=\"\" part=\"0\" name=\"\" file=\""
                    << seqName_(i) << ".vtp" << "\"/> \n";
        }

        pvdFile << "</Collection> \n"
                << "</VTKFile> \n" << std::flush;
    }

    //! add cell data vector
    void addCellData(const std::vector<double>& cellData, const std::string& name)
    { vtkWriter_.addCellData(cellData, name); }

private:
    std::string seqName_(std::size_t count) const
    {
        std::stringstream n;
        n << name_ << "-" << std::setfill('0') << std::setw(5) << count;
        return n.str();
    }

    RBCNetworkVTKWriter<RBCNetwork> vtkWriter_;
    const std::string name_;
    std::vector<double> timeSteps_;
};

} // end namespace Dumux

#endif
