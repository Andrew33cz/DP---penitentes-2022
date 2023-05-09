/* ------------------------------------
*      bc. Ondrej Cech (xcecho06), 2022
*  ------------------------------------
*/

#include "functions.hpp"

using namespace std;



int main(int argc, char *argv[])
{
    if(argc!=5){
        cerr<<"Expected parameters: input file name, length of one step (seconds), number of steps, number of rays."<<endl<<"For example penitentes.exe ./Test/Test.ply 3600 24 1000 "<<endl;
        return 1;
    }
    else if(atoi(argv[2])<1){
        cerr<<"Argument 2 (step length) must be natural number (greater than 0)."<<endl;
        return 1;
    }
        else if(atoi(argv[3])<1){
        cerr<<"Argument 3 (number of steps) must be natural number (greater than 0)."<<endl;
        return 1;
    }
        else if(atoi(argv[4])<1){
        cerr<<"Argument 4 (number of rays) must be natural number (greater than 0)."<<endl;
        return 1;
    }
    string name(argv[1]);
    cout<<name<<" "<<atoi(argv[2])<<" "<<atoi(argv[3])<<" "<<atoi(argv[4])<<" "<<endl;
    simulate(name,atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
    //smtest(name,atoi(argv[4]));
    //cout<<derviaceH({1.0,0.0,-0.5},{562.0,25.2,35.8},{-1.0,0.0,1.0})<<endl;
    return 0;
}
