#include <iostream>
#include <armadillo>
#include<fstream>
#include <string>

using namespace std;
//using namespace arma;



void get_noise_cov(string filename)
{
    filename = string(filename);
    ifstream infile;
    infile.open (filename);
    for( string line; getline( infile, line ); )
    {
      cout << line[65:75] << endl;
    }
    //return N_cov;

}

int main()
{
    get_noise_cov("90GHz.txt");
    return 0;
}

