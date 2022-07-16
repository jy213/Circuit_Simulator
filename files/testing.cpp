#include <iostream>
#include <vector>
#include <cmath>
#include<fstream>
#include<cstdlib>
#include<sstream>
using namespace std;

void multiplier(double &value,char m){
        if(m == 'p'){
            value *= pow(10, -12);
        }
        if(m == 'n'){
            value *= pow(10, -9);
        }
        if(m == 'u'){
            value *= pow(10, -6);
        }
        if(m == 'm'){
            value *= pow(10, -3);
        }
        if(m == 'k'){
            value *= pow(10, 3);
        }
        if(m == 'g'){
            value *= pow(10, 6);
        }    
        if(m == 'G'){
            value *= pow(10, 9);
        }
    }

string getword(string tmp, int wd){ //the leftmost one of the line is word 1, not word 0
    string str[6];
    istringstream is(tmp);
    for(int i=0; i<6; i++){
        is >> str[i];
    }
    return str[wd-1];
}

string lastword(string tmp){
    string lwd;
    for(int i=tmp.size()-1; tmp[i]!=' '; i--){
        lwd += tmp[i];
    }
    string s(lwd.rbegin(),lwd.rend());
    return s;
}

double ConvertFromString(string str) { //string to double
    istringstream iss(str);
    double m;
    if (iss >> m)
    return m;
    return 0.0;
}

void mul(double &val, string s){
    if(s[s.size()-1] == 'g'){
        string value;
        for(int j=0; j<s.size() - 3; j++){
            value.push_back(s[j]);
        }
        val= ConvertFromString(value);
        multiplier(val, 'g');
    }
    else if(s[s.size()-1] == 'p' || s[s.size()-1] == 'n' || s[s.size()-1] == 'u' || s[s.size()-1] == 'm'
        || s[s.size()-1] == 'k' || s[s.size()-1] == 'G'){
        string value;
        for(int j=0; j<s.size(); j++){
            value.push_back(s[j]);
        }
        val= ConvertFromString(value);
        multiplier(val, s[s.size()-1]);
    }   
    else{
        val = ConvertFromString(s);
    }
}

vector<double> ac (vector<string> input){
    vector<double> ac;
    string p = getword(input[input.size()-2], 3);
    double points;
    mul(points, p);
    ac.push_back(points);
    string startf = getword(input[input.size()-2], 4);
    double start;
    mul(start, startf);
    ac.push_back(start);
    string stopf = getword(input[input.size()-2], 5);
    double stop;
    mul(stop, stopf);
    ac.push_back(stop);
    return ac;
}

std::vector<double> frequency (std::vector<double> ac){
    std::vector<double> tmp;
    double n = log10(ac[2] / ac[1]);
    double fs = ac[1];
    double fn = 0;
    for(int i = 0; i <= n; i++){
        for(int g = 0; g < ac[0]; g++ ){
           fn = fs * pow(10, g/ac[0]);
           tmp.push_back(fn);
        }
        fs = 10 * fs;   
    }
    return tmp;
}
int main(){
    vector<double> ac;
    vector<double> fre;
    ac.push_back(1);
    ac.push_back(10);
    ac.push_back(100000);
    fre = frequency(ac);
    for(int i=0;i<fre.size();i++){
        cout<<fre[i]<<endl;
    }
}