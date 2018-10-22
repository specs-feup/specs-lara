#include <iostream>
#include <H5Cpp.h>

#include "Data.h"

#include "CompType.h"

#define DATASET "TheDataset"
#define H5FILE "data.h5"

using namespace H5;
using namespace H5EX;

void writeData(const Data &data)
{
      Data dataArr[] = {data};
      
      /*
       * Create the data space.
       */
      hsize_t dim[] = {1};   /* Dataspace dimensions */
      DataSpace space( 1, dim );
      /*
       * Create the file.
       */
      H5File* file = new H5File( H5FILE, H5F_ACC_TRUNC );
       
      CompType type = hdf5type::DataType::GetCompType(); // Obtain compound type mapping

      /*
       * Create the dataset.
       */
      DataSet* dataset;
      dataset = new DataSet(file->createDataSet(DATASET, type, space));
      /*
       * Write data to the dataset;
       */
      dataset->write( dataArr, type);
      /*
       * Release resources
       */
      delete dataset;
      delete file;

}


Data readData()
{
      /*
       * Open the file and the dataset.
       */
      H5File* file = new H5File(H5FILE , H5F_ACC_RDONLY );
      DataSet* dataset = new DataSet (file->openDataSet( DATASET ));
      
      /*
       * Obtain a compound datatype
       */
      CompType type = hdf5type::DataType::GetCompType();
      
      /*
       * Read a dataset
       */
      Data dataArr[1];
      dataset->read( dataArr, type );
      
      delete dataset;
      delete file;

      return dataArr[0];
}


int main()
{
   
  // Create a structure instance and store it
  Data d;
  d.id = 1000;
  d.lat = 14.345;
  d.lon = 48.235;
  d.cls = 'a';
  d.velocity = 234.1;

  // Write the structure
  writeData(d);


  // Read it back
  Data dr = readData();


  // Print its content
  std::cout << "Got Data: ";
  std::cout << "id: " << dr.id;
  std::cout << " Coords: " << dr.lat << "," << dr.lon;
  std::cout << " Class: " << dr.cls;
  std::cout << " Velocity: " << dr.velocity << std::endl;

  return 0;
}
