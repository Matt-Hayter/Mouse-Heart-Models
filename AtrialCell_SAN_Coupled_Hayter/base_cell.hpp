#ifndef PARENT_CELL_HPP
#define PARENT_CELL_HPP

#include <fstream>

//Abstract base class for all cell models. All pure virtual functions mmust be implemented by each cell type

class base_cell
{
    public:
        virtual ~base_cell() {}

        virtual void output_config(std::ofstream &) = 0;
        virtual void create_data_file() = 0;
        virtual void create_measurement_file() = 0;
        virtual void set_initial_states() = 0;

        virtual void store_variables() = 0;
        virtual void ODEs(const double &) = 0;
        virtual void Euler_method() = 0;
        virtual void measurements(const double &, int &) = 0;

        virtual void output_data(const double &) = 0;
        virtual void output_measurements(const int &) = 0;
        virtual void output_final_states() = 0;

        virtual const double &get_potential() = 0;
        virtual void set_I_j(const double &) = 0;
};

#endif