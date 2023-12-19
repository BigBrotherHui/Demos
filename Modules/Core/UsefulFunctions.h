#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <Eigen/Eigen>
#include <vtkMatrix4x4.h>
#include <vtkSmartPointer.h>
using namespace std;
using namespace Eigen;
//注意：通过Map构造的矩阵相当于引用，会改变源数据的值，是浅拷贝的一种实现形式
MatrixXd vector2MatrixXd(vector<double> v,int row)
{
    return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(v.data(), row, v.size() / row);
}
void Eigen2VtkMatrix(const Eigen::Matrix4d &mt,vtkSmartPointer<vtkMatrix4x4> vmt)
{
    memcpy(vmt->GetData(), mt.data(), sizeof(double) * 16);
    vmt->Transpose();
}
//shallow copy internal data
void VtkMatrix2Eigen(vtkSmartPointer<vtkMatrix4x4> vmt,Eigen::Matrix4d &mt)
{
    mt = Map<Matrix4d>(vmt->GetData()).transpose();
    /*vtkSmartPointer<vtkMatrix4x4> vmt = vtkSmartPointer<vtkMatrix4x4>::New();
	vmt->SetElement(1, 3, 24);
	Eigen::Matrix4d mt;
	VtkMatrix2Eigen(vmt, mt);
	std::cout << mt;*/
}
MatrixXd readCsv(string fileToOpen)
{
    vector<double> matrixEntries;
    // in this object we store the data from the matrix
    ifstream matrixDataFile(fileToOpen);

    // this variable is used to store the row of the matrix that contains commas 
    string matrixRowString;

    // this variable is used to store the matrix entry;
    string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;

    getline(matrixDataFile, matrixRowString);
    matrixRowString = "";

    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

        while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
        {
            // matrixEntries.push_back(matrixEntries.cast<double>(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
            matrixEntries.push_back(stod(matrixEntry));
        }
        matrixRowNumber++; //update the column numbers
    }

    // here we convert the vector variable into the matrix and return the resulting object, 
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
    return vector2MatrixXd(matrixEntries, matrixRowNumber);

}

void readCsv(const string& name, vector<double>& X, vector<double>& Y)
{
    ifstream csv_data(name, ios::in);
    if (!csv_data.is_open())
    {
        cout << "Error: opening file fail" << std::endl;
        exit(1);
    }

    string line;
    vector<string> lines;
    while (getline(csv_data, line))
        lines.emplace_back(line);

    int number = lines.size();
    X.resize(number);
    Y.resize(number);
    for (int i = 0; i < number; i++)
    {
        istringstream iss(lines[i]);
        string token;

        getline(iss, token, ',');
        X[i] = atof(token.c_str());

        getline(iss, token, ',');
        Y[i] = atof(token.c_str());
    }

    csv_data.close();
}