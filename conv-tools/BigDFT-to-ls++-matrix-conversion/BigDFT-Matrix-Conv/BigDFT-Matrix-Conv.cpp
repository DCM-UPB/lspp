//Programm for converting BigDFT Matrixformat to linear Scaling CRS-format.
//Author: H. Wiebeler

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;
using std::vector;

  int main(int argc, char *argv[]){
    if (argc < 2){
      cout << "Syntax is: BD-Matriz-Read BigDFT Matriz file!" << endl;
      return 1;
      }
    
    ifstream DataIn;
    stringstream data;
    ofstream DataOut;

    string line;
    string answer;
    int irow;     
    int irow_old = 0; 
    int jcolum;     
    int jmax = 0;  //max value for j
    int valzahl = 0;  //number of values
    int anzahl_rp = 0;   //number of rowpointers
    vector<int> rowpointers(31623); //assumption sparsity < 10 % 
    vector<int> indices(100000000); //therefore 100.000.000 max  
    vector<double> values(100000000);
    int rowpointer = 0;
    double value;

    DataIn.open (argv[1],ios::in); 
    cout << "Reading matrix..." << endl;
    while (getline(DataIn, line)) {
      
     
      if (line[0] == '#')
        continue;
      
      data << line; 
      data >> irow;
      data >> jcolum;
      data >> value;            //Read and saving of BigDFT-Matrix 
      data.str("");

      if (jcolum > jmax)
        jmax = jcolum;

      if (value == 0)
        continue;
      if (irow_old < irow){
        rowpointers[anzahl_rp] = rowpointer;
        anzahl_rp++;
      
      if (irow_old +1 != irow){
        cout << "Error, row is missing, please change the code" << endl;
        return 11;}  
       }
      rowpointer++;
      indices[valzahl] = jcolum - 1;          //take care of Matrix-shift BigDFT 1,1... lin_Scaling 0,0...!
      values[valzahl] = value;
      valzahl++;
      irow_old = irow;

      if (anzahl_rp == 31623){
        cout << "Error, matrix is to big, please change the code" << endl;     //Error in case of to big matrices
        return 1;}  
      continue; 
        }
    rowpointers[anzahl_rp] = rowpointer;
 

   if (irow < jmax) {
    cout << "Warning the last or more rows seems to have no values!!!" << endl; 
    cout << "Proceed (y/n)?" << endl;      //Error if rows are missing.
    cin >> answer;

    if (answer == "y")
      irow = jmax;                        //Makes sure, that irow has the same dimension as the quad. matrix.
          
    
    else
     return 2;
   }
   
    
    cout << "Writing Matrix." << endl;
    DataOut.open ("Matrix",ios::out);     //Matrix Output
    DataOut << "rows " << irow << endl;
    DataOut << "cols " << irow << endl;

    DataOut << "rowpointers " << irow + 1 << endl;
    for (int m = 0; m <= irow ; m++)              //rowpointers
      DataOut << rowpointers[m] << " ";
    DataOut << endl;
    
    DataOut << "indices " << valzahl << endl;
    for (int n = 0; n < valzahl; n++)        //Indices
      DataOut << indices[n] << " ";
    DataOut << endl;

    DataOut << "values " << valzahl << endl;
    for (int o = 0; o < valzahl; o++)
      DataOut << setprecision(12) << scientific << values[o] << " " ; //Values
    DataOut << endl;

          
    return 0;
}
