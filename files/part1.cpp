#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>
#include<fstream>
#include<cstdlib>
#include<sstream>
#include<math.h>
#define PI 3.14159265 
const double VT = 0.025;
//bjt
// pick 2N2907 PNP VAF = 120 BF = 250 Is = 1E-14 BR = 3;
// pick 2N2222 NPN VAF = 100 BF = 200 Is = 1E-14 BR = 3;
double VAnpn = 120;
double VApnp = 100;
double BFnpn = 200;
double BFpnp = 250;
double BR = 3;
double Is = pow(10, -14);
double Isc = ((1+BR) / BR) * pow(10, -14);
//mosfet
// for NMOS vt = 2.9 kp = 0.9m l = 10u w = 100u lambda = 0.01 VA = 100;
// for PMOS vt = -2.5 kp = 0.4m l = 10u w = 100u lambda = 0.03 VA = 300;
double Kn = 0.009;
double Kp = 0.004;
double vtn = 2.9;
double vtp = -2.5;
double VAn = 100;
double VAp = 300;
class general{
public: 
    virtual std::complex<double> conductance (double omega) {return 0;}
    virtual int nodes() const = 0;
    general(int n1, int n2, std::string t) : node_1(n1), node_2(n2), type(t){}
    std::vector<int> get_node(){
        std::vector<int> v;
        if(node_1 > node_2){
            v.push_back(node_2);
            v.push_back(node_1);
        }
        else{
            v.push_back(node_1);
            v.push_back(node_2);
        }
        return v;
    }
    virtual std::string get_type (void){
        return type;
    }
    virtual std::complex<double> get_volt(void){
        return 0;
    }
    virtual std::complex<double> get_current(void){
        return 0;
    }
    virtual std::vector<int> get_polarity(void){
        std::vector<int> g;
        g.push_back(node_1);
        g.push_back(node_2);
        return g;
    }
    virtual void change_node(int n1, int n2){
        node_1 = n1;
        node_2 = n2;
    }
    virtual~ general(){ };
protected:
    int node_1;
    int node_2;
    std::string type;
};
class Resistor : public general{
public:
    Resistor(int n1, int n2, double value, std::string t) : general(n1, n2, t), resistance(value){} 
    std::complex<double> conductance(double omega) {
        return (1 / resistance);
    }
    int nodes()const{
        return 2;
    }
    double real_cond(){
        return (1/resistance);
    }
private:
    double resistance;
};
class Capacitor : public general{
public:
    Capacitor(double c, int n1, int n2, std::string t) : capacitance(c) , general(n1, n2, t){ }  
    std::complex<double> conductance(double omega) {
        std::complex<double> tmp;
        tmp.real(0);
        tmp.imag(omega * capacitance);
        return tmp;
    }
    int nodes()const {
        return 2;
    }
private:
    double capacitance;
};
class Inductor : public general{
public:
    Inductor(double l, int n1, int n2, std::string t) : general(n1,n2, t), inductance(l) {     }  
    std::complex<double> conductance(double omega){
        std::complex<double> tmp;
        tmp.real(0);
        tmp.imag(- 1 / (omega * inductance));
        return tmp;
    }
    int nodes() const {
        return 2;
    }
private:
    double inductance;
};
class Diode : public general { //model 1N914
public:
    Diode(int n1, int n2, std::string t) : general(n1, n2, t){ }
    std::complex<double> conductance(double omega)const{
        return 0;
    }
    double id (double v1, double v2){
        double tmp = 0;
        tmp = pow(10, -12) * (exp((v1 - v2)/VT)-1);
        return tmp;
    }
    int nodes() const{
        return 2;
    }
    double g (double v1, double v2){
        double tmp = pow(10, -12) * exp((v1 - v2)/VT) / VT;
        return tmp;
    }
    std::vector<int> get_polarity(){
        std::vector<int> tmp;
        tmp.push_back(node_1);
        tmp.push_back(node_2);
        return tmp;
    }
private:
};
class voltsrc : public general{
public:
    voltsrc(double v, int n1, int n2, std::string t, std::string n) : voltage(v), general(n1, n2, t), amp(0), phase(0), name(n) {}
    voltsrc(double a, double p, int n1, int n2, std::string t, std::string n) : voltage(0), amp(a), phase(p), general(n1, n2, t), name(n){}
    std::complex<double> conductance(double omega)const{
        return 0;
    }
    std::complex<double> get_volt () {
        if (amp == 0 && voltage != 0){
            return voltage;
        }
        else {
            std::complex<double> tmp;
            tmp.real( amp * cos(phase * PI / 180.0));
            tmp.imag(amp * sin(phase * PI /180.0));
            return tmp;
        }
    }
    std::string dc_ac() {
        if (amp == 0 && voltage != 0){
            return "dc";
        }
        if (amp != 0 && voltage == 0){
            return "ac";
        }
        else{
            return 0;
        }
    }
    int nodes() const{
        return 2;
    }
    std::vector<int> get_polarity(){
        std::vector<int> g;
        g.push_back(node_1);
        g.push_back(node_2);
        return g;
    }
    std::string search_name(){
        return name;
    }
private:
    double voltage;
    double amp;
    double phase;
    std::string name;
};
class currsrc : public general{
public:
    currsrc(double i, int n1, int n2, std::string t, std::string n) : amp(0), phase(0), current(i) , general(n2 , n1, t), name(n){}
    currsrc(double a, double p, int n1, int n2, std::string t, std::string n) : current(0), amp (a) , phase(p), general(n2, n1, t) ,name(n) {         }
    std::complex<double> conductance(double omega)const{
        return 0;
    }
    std::complex<double> get_current (){
        if (amp == 0 && current != 0){
            return current;
        }
        else {
            std::complex<double> tmp;
            tmp.real(amp * cos(phase * PI / 180.0));
            tmp.imag(amp * sin(phase * PI /180.0));
            return tmp;
        }
    }
    std::string dc_ac(){
        if (amp == 0 && current != 0){
            return "dc";
        }
        if (amp != 0 && current == 0){
            return "ac";
        }
        else{
            return 0;
        }
    }
    int nodes() const{
        return 2;
    }
    std::vector<int> get_polarity(){
        std::vector<int> g;
        g.push_back(node_1);
        g.push_back(node_2);    
        return g;
    }
    std::string search_name(){
        return name;
    }
private:
    double current;
    double amp;
    double phase;
    std::string name;
};
class bjt: public general{
public:
    bjt(int n1, int n2, int n3, std::string tmp, std::string t) : general(n1, n2, t), type_b(tmp){
        collector = n1;
        base = n2;
        emitter = n3;
    }
    std::complex<double> conductance(double omega)const{
            return 0;
    }
    void initialize_bjt (double v1 , double v2, double v3){
        vb = v1;
        ve = v2;
        vc = v3;
        if(type_b == "NPN"){
            if(vb>ve){
                if(vc>vb){
                    state = "A";
                    Ic = Is * exp((vb-ve)/VT) * (1+ (vc-ve)/VAnpn);
                    Ib = Ic / BFnpn;
                    Ie = (BFnpn + 1)/BFnpn * Ic;
                }
                else{
                    state = "S";
                    Ic = Is * (exp((vb-ve)/VT)-exp((vb-vc)/VT)) - Isc * (exp((vb-vc)/VT)-1);
                    Ib = Is * (exp((vb-ve)/VT)-1)/BFnpn + Isc * (exp((vb-vc)/VT)-1);
                    Ie = Is * (exp((vb-ve)/VT)-exp((vb-vc)/VT)) + Is * (exp((vb-ve)/VT)-1)/BFnpn;
                }
            }
            else{
                if(vc>vb){
                    state = "C";
                }
                else{
                    state = "R";
                }
                Ic = Is * (exp((vb-ve)/VT)-exp((vb-vc)/VT)) - Isc * (exp((vb-vc)/VT)-1);
                Ib = Is * (exp((vb-ve)/VT)-1)/BFnpn + Isc * (exp((vb-vc)/VT)-1);
                Ie = Is * (exp((vb-ve)/VT)-exp((vb-vc)/VT)) + Is * (exp((vb-ve)/VT)-1)/BFnpn;
            }
        }
        if(type_b == "PNP"){
            if(vb>ve){
                if(vc>vb){
                    state = "R";
                }
                else{
                    state = "C";
                    Ic = 0;
                    Ib = 0;
                }
                Ic = Is * (exp((ve-vb)/VT)-exp((vc-vb)/VT)) - Isc * (exp((vc-vb)/VT)-1);
                Ib = Is * (exp((ve-vb)/VT)-1)/BFpnp + Isc * (exp((vc-vb)/VT)-1);
                Ie = Is * (exp((ve-vb)/VT)-exp((vc-vb)/VT)) + Is * (exp((ve-vb)/VT)-1)/BFpnp;
            }
            else{
                if(vc>vb){
                    state = "S";
                    Ic = Is * (exp((ve-vb)/VT)-exp((vc-vb)/VT)) - Isc * (exp((vc-vb)/VT)-1);
                    Ib = Is * (exp((ve-vb)/VT)-1)/BFpnp + Isc * (exp((vc-vb)/VT)-1);
                    Ie = Is * (exp((ve-vb)/VT)-exp((vc-vb)/VT)) + Is * (exp((ve-vb)/VT)-1)/BFpnp;
                }
                else{
                    state = "A";
                    Ic = Is * exp((ve-vb)/VT) * (1+ (ve-vc)/VApnp);
                    Ib = Ic / BFpnp;
                    Ie = (BFpnp + 1)/BFpnp * Ic;
                }
            }
        }
    }
    std::vector<int> get_polarity (){
        std::vector<int> tmp;
        tmp.push_back(collector);
        tmp.push_back(base);
        tmp.push_back(emitter);
        return tmp;
    }
    std::complex<double> rbe (){
        if(type_b == "NPN"){
            return (VT * BFnpn) / Ic;
        }
        else{
            return (VT * BFpnp) / Ic;
        }
    }
    std::complex<double> r0 (){
        if(type_b == "NPN"){
            return VAnpn / Ic;
        }
        else{
            return VApnp / Ic;
        }
    }
    std::complex<double> gm (){
        return Ic / VT;
    }
    int nodes() const{
        return 3;
    }
    double b_b(){
        if(state != "A"){
            if(type_b == "NPN"){
                return Is * exp((vb-ve)/VT) / (VT*BFnpn) + Isc * exp((vb-vc)/VT) / VT;
            }
            else{
                return -Is * exp((ve-vb)/VT) / (VT*BFpnp) - Isc * exp((vc-vb)/VT) / VT;
            }
        }
        else{
            if(type_b =="NPN"){
                return Ib / VT;
            }
            else{
                return -Ib / VT;
            }
        }
    }
    double b_c(){
        if(state != "A"){
            if(type_b == "NPN"){
                return -Isc * exp((vb-vc) / VT) / VT;
            }
            else{
                return Isc * exp((vc-vb) / VT) / VT;
            }
        }
        else{
            if(type_b == "NPN"){
                return Is * exp((vb-ve) / VT) / (VAnpn*BFnpn);
            }
            else{
                return -Is * exp((ve-vb) / VT) / (VApnp*BFpnp);
            }
        }
    }
    double b_e(){
        if(state != "A"){
            if(type_b == "NPN"){
                return -Is * exp((vb-ve) / VT) / (VT * BFnpn);
             }
            else{
                return Is * exp((ve-vb) / VT) / (VT*BFpnp);
            }
        }
        else{
            if(type_b == "NPN"){
                return -(Ib / VT + Is * exp((vb-ve)/VT) / (BFnpn*VAnpn));
            }
            else{
                return Ib / VT + Is * exp((ve-vb)/VT) / (BFpnp*VApnp);
            }
        } 
    }
    double c_c(){
        if(state == "A"){
            if(type_b == "NPN"){
                return Is * exp((vb-ve)/VT) / VAnpn;
            }
            else{
                return - Is * exp((ve-vb)/VT)/ VApnp;
            }  
        }
        else{
            if(type_b == "NPN"){
                return ((Is + Isc)* exp((vb-vc)/VT))/VT;
            }
            else{
                return - ((Is + Isc)* exp((vc-vb)/VT))/VT;
            }
        }
    }
    double c_e(){
        if(state == "A"){
            if(type_b == "NPN"){
                double x =0;
                x = - Is * exp((vb-ve)/VT) / VAnpn - Ic/VT;
                return x;
            }
            else{
                double y =0;
                y = Is * exp((ve-vb)/VT) / VApnp + Ic/VT;
                return y;
            }  
        }
        else{
            if(type_b == "NPN"){
            return -Is* exp((vb-ve)/VT)/VT;
            }
            else{
            return Is * exp((ve-vb)/VT)/VT;
            }
        }
    }
    double c_b(){
        if(state == "A"){
            if(type_b == "NPN"){
                return Ic/VT;
            }
            else{
                return - Ic/VT;
            }
        }
        else{
            if(type_b == "NPN"){
                return (Ic-Isc)/VT;
            }
            else{
                return (Isc-Ic)/VT;
            }
        }
    }
    double e_e(){
        if(state == "A"){
            if(type_b == "NPN"){
                double x1 =0;
                x1 = - Is * exp((vb-vc)/VT) / VAnpn - Ic/VT;
                x1 = x1 * (BFnpn+1)/BFnpn;
                return x1;
            }
            else{
                double y1 =0;
                y1 = Is * exp((vc-vb)/VT) / VApnp + Ic/VT;
                y1 = y1 * (BFpnp + 1)/ BFpnp;
                return y1;
            }    
        }
        else{
            if(type_b == "NPN"){
                return -(Is + Is/BFnpn) * exp((vb-ve)/VT)/VT;
            }
            else{
                return (Is/BFpnp+Is) * exp((ve-vb)/VT)/VT;
            }
        }
    }
    double e_b(){
        if(state == "A"){
            if(type_b == "NPN"){
                return Ie/VT;
            }
            else{
            return -Ie/VT;
            }
        }
        else{
            if(type_b == "NPN"){
                return (Ie + Is/BFnpn)/VT;
            }
            else{
                return -(Ie + Is/BFpnp)/VT;
            }  
        }
    }
    double e_c(){
        if(state == "A"){
            if(type_b == "NPN"){
                return Is * exp((vb-ve)/VT) / VAnpn * (BFnpn+1)/BFnpn;
            }
            else{
                return -Is * exp((ve-vb)/VT)/ VApnp * (BFpnp+1)/BFpnp;
            }  
        }
        else{
            if(type_b == "NPN"){
                return Is * exp((vb-vc)/VT)/VT;
            }
            else{
                return -Is * exp((vc-vb)/VT)/VT;
            }
        }
    }
    double ib(){
        return Ib;
    }
    double ic(){
        return Ic;
    }
    double ie(){
        return Ie;
    }
    int b(){
        return base;
    }
    int e(){
        return emitter;
    }
    int c(){
        return collector;
    }
    void change_node( int n1, int n2, int n3){
        collector = n1;
        base = n2;
        emitter = n3;
    }
    std::string getmodel(){
        return type_b;
    }
private:
    int base;
    int emitter;
    int collector;
    std::string type_b;
    double vb;
    double ve;
    double vc;
    double Ic;
    double Ib;
    double Ie;
    std::string state;
};
class mosfet : public general{
public:
    mosfet(int n1, int n2, int n3, std::string tmp, std::string t) : general(n1 , n2, t), type(tmp){
        drain = n1;
        gate = n2;
        source = n3;
    }
    std::string getmodel(){
        return type;
    }
    double r0(){
        if(type == "NMOS"){
            return  VAn/ Id;
        }
        else{
            return  VAp / Id;
        }
    }
    double gm(){
        if(type == "NMOS"){
            return 2 * sqrt(Kn * Id);
        }           
        else{
            return 2 * sqrt(Kp * Id);
        }
    }
    std::complex<double> conductance(double omega)const{
        return 0;
    }
    int nodes() const{
        return 3;
    }
    std::vector<int> get_polarity (){
        std::vector<int> tmp;
        tmp.push_back(gate);
        tmp.push_back(source);
        tmp.push_back(drain);
        return tmp;
    }
    void initialize_mos (double v1 , double v2, double v3){
        vg = v1;
        vs = v2;
        vd = v3;
        double tmp_vt = 0;
        if(type == "NMOS"){ 
            if(vg-vs >= vtn){
                if(vd >= vg-vtn){
                    state = "S";
                    Id = Kn * pow((vg-vs-vtn), 2) * (1 + (vd- vs)/VAn);
                }
                else{
                    state = "T";
                    Id = Kn * (2 * (vg-vs-vtn) * (vd-vs) - pow(vd-vs, 2));
                }
            }
            else{
                state = "OFF";
                Id = 0;
            }
        }
        else{
            if(vg-vs < vtp){
                if(vd <= vg-vtp){
                    state = "S";
                    Id = Kp * pow((vg-vs-vtp), 2) * (1+ (vs - vd)/VAp);
                }   
                else{
                    state = "T";
                    Id = Kp * (2 * (vg-vs-vtp) * (vd-vs) - pow(vd-vs, 2));
                }
            }
            else{
                state = "OFF";
                Id = 0;
            }
        }
    }
    double Id_d (){
        // provided the direction of the current    
        if(state == "T"){
            if(type == "NMOS"){
                return 2*Kn*(vg-vtn-vd);
            }
            else{
                return 2*Kp*(vg-vd+abs(vtp));
            }
        }
        else{
            if(type == "NMOS"){
                 return Kn * pow((vg-vs-vtn), 2) / VAn;
            }
            else{
                return -Kp * pow((vs-vg-abs(vtp)), 2) / VAp;
            }
        }
    }
    double Id_s(){
        if(state == "T"){
            if(type == "NMOS"){
                return 2 * Kn * (vs - vg + vtn);
            }
            else{
                return 2 * Kp * (vs - vg - abs(vtp));
            }
        }
        else{
             if(type == "NMOS"){
                return -2*Kn*(vg-vs-vtn)*(1+(vd-vs)/VAn)+Kn*pow((vg-vs-vtn), 2)*(-1/VAn);
            }
            else{
                return 2*Kp*(vs-vg-abs(vtp))*(1+(vs-vd)/VAp)+Kp*pow((vs-vg-abs(vtp)), 2)*(1/VAp);
            }
        }
    }
    double Id_g(){
        if(state == "T"){
            if(type == "NMOS"){
                return 2 * Kn * (vd - vs);
            }
            else{
                return -2 * Kp * (vs - vd);
            }
        }
        else{
            if(type == "NMOS"){
                return 2 * Kn * (vg - vs -vtn) * (1 + (vd - vs)/VAn);
            }
            else{
                return -2 * Kp * (vs - vg - abs(vtp)) * (1 + (vs - vd)/VAp);
            }
        }
    }
    double id(){
        return Id;
    }
    double vt(){
        if(type == "NMOS"){
            return vtn;
        }
        else{
            return vtp;
        }
    }
    int g(){
        return gate;
    }
    int s(){
        return source;
    }
    int d(){
        return drain;
    }
    void change_node( int n1, int n2, int n3){
        drain = n1;
        gate = n2;
        source = n3;
    }
private:
    int gate;
    int source;
    int drain;
    double vg;
    double vs;
    double vd;
    double Id;
    std::string type;
    std::string state;
};
class v_currsrc : public general{
public:
    v_currsrc(int n1, int n2, int n3, int n4, double t, std::string tm) : general(n1, n2, tm), control_p(n3), control_n(n4), trans(t){}
    int nodes() const{
        return 4;
    }
    std::complex<double> get_i (std::complex<double> v1, std::complex<double> v2){ // v1 = v+*
        return (v1 - v2) * trans;
    }
    std::complex<double> conductance(double omega){
        return trans;
    }
    std::vector<int> get_polarity(){
        std::vector<int> tmp;
        tmp.push_back(node_1);
        tmp.push_back(node_2);
        tmp.push_back(control_p);
        tmp.push_back(control_n);
        return tmp;
    }
private:
    int control_p;
    int control_n;
    double trans;
};

std::vector <std::string> ReadInput (std::string filename){ //read file, store in vector of type string
    std::ifstream infile;
    infile.open(filename);
    if(!infile.is_open()){
        std::cout << "error opening file" << std::endl;
    }
    std::vector <std::string> input;
    std::string tmp;
    while(getline(infile, tmp)){
        input.push_back(tmp);
    }
    return input;
    infile.close();
}

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
std::string getword(std::string tmp, int wd){ //the leftmost one of the line is word 1, not word 0
    std::string str[6];
    std::istringstream is(tmp);
    for(int i=0; i<6; i++){
        is >> str[i];
    }
    return str[wd-1];
}
std::string lastword(std::string tmp){
    std::string lwd;
    for(int i=tmp.size()-1; tmp[i]!=' '; i--){
        lwd += tmp[i];
    }
    std::string s(lwd.rbegin(),lwd.rend());
    return s;
}
double ConvertFromString(std::string str) { //string to double
    std::istringstream iss(str);
    double m;
    if (iss >> m)
    return m;
    return 0.0;
}
void multi(double &value, std::string s){
//for phase calculation
    if(s[s.size()-1] == ')'){
        char seclast = s[s.size()-2];
        bool sec = false;
        if(seclast == 'p' || seclast == 'n' || seclast == 'u' || seclast == 'm' || seclast == 'k' || seclast == 'g' || seclast == 'G'){
            sec = true; //phase is followed by a multiplier
        }
        if(sec == true){
            if(seclast == 'g'){
                value = ConvertFromString(s.substr(0, s.size()-4));
                multiplier(value, seclast);
            }
            else{
                value = ConvertFromString(s.substr(0, s.size()-2));
                multiplier(value, seclast);
            }
        }
        else{
            value = ConvertFromString(s.substr(0, s.size()-1));
        }
    }
//for other calculation (including ac and dc source value)
    else{
        bool m = false;
        char last = s[s.size()-1];
        //see if there is a multiplier
        if(last == 'p' || last == 'n' || last == 'u' || last == 'm' || last == 'k' || last == 'g' || last == 'G'){
            m = true;
        }
        //consider AC case
        if(s[0] == 'A' && m == true){ //AC input with multiplier
            if(last == 'g'){
                value = ConvertFromString(s.substr(3, s.size()-6));
                multiplier(value, last);
            }
            else{
                value = ConvertFromString(s.substr(3, s.size()-4));
                multiplier(value, last);
            }
        }
        if(s[0] == 'A' && m != true){ //AC input without multiplier
            value = ConvertFromString(s.substr(3));
        }
        //consider DC case
        if(s[0] != 'A' && m == true){ //DC input with multiplier
            if(last == 'g'){
                value = ConvertFromString(s.substr(0, s.size()-3));
                multiplier(value, last);
            }
            else{
                value = ConvertFromString(s.substr(0, s.size()-1));
                multiplier(value, last);
            }
        }
        if(s[0] != 'A' && m != true){ //DC input without multiplier
            value = ConvertFromString(s.substr(0));
        }
    }
}
std::vector<double> ac (std::vector<std::string> input){
    std::vector<double> ac;
    std::string p = getword(input[input.size()-2], 3);
    double points=0;
    multi(points, p);
    ac.push_back(points);
    std::string startf = getword(input[input.size()-2], 4);
    double start=0;
    multi(start, startf);
    ac.push_back(start);
    std::string stopf = getword(input[input.size()-2], 5);
    double stop=0;
    multi(stop, stopf);
    ac.push_back(stop);
    return ac;
}

void setting (std::vector <std::string> input, std::vector<general*> &component, int &max_node){ //split lines into parts, store in vector of type component
    general* tmp_a;
    for(int i=0; i<input.size()-2; i++){
        if(input[i][0] != '*'){
            double value=0;
            std::string node1, node2;
            int n1, n2;
            node1 = getword(input[i], 2);
            node2 = getword(input[i], 3);
            if(node1 == "0"){
                n1 = 0;
            }
            else{
                n1 = std::stoi(node1.substr(1));
            }
            if(n1 > max_node){
                max_node = n1;
            }
            if(node2 == "0"){
                n2 = 0;
            }
            else{
                n2 = std::stoi(node2.substr(1));
            }
            if(n2 > max_node){
                max_node = n2;
            }
            if(getword(input[i], 1)[0] == 'V'){
                if(getword(input[i], 4)[0] == 'A'){ //value or function
                    multi(value, getword(input[i], 4));
                    double phase=0;
                    multi(phase, getword(input[i], 5));
                    tmp_a = new voltsrc (value, phase, n1 , n2, "V", getword(input[i] , 1));
                    component.push_back(tmp_a);
                } 
                else{
                    multi(value, getword(input[i], 4));
                    tmp_a = new voltsrc(value, n1, n2, "V", getword(input[i], 1));
                    component.push_back(tmp_a);
                }
            }
            if(getword(input[i], 1)[0] == 'I'){
                if(getword(input[i], 4)[0] == 'A'){ //value or function
                    multi(value, getword(input[i], 4));
                    double phase=0;
                    multi(phase, getword(input[i], 5));
                    tmp_a = new currsrc (value, n1, n2, "I", getword(input[i], 1));
                    component.push_back(tmp_a);
                }             
                else{
                    multi(value, getword(input[i], 4));
                    tmp_a = new currsrc(value, n1, n2, "I", getword(input[i], 1));
                    component.push_back(tmp_a);
                } 
            }
            if(getword(input[i], 1)[0] == 'R'){
                multi(value, getword(input[i], 4));
                tmp_a = new Resistor (n1, n2, value, "R");
                component.push_back(tmp_a);
            }
            if(getword(input[i], 1)[0] == 'C'){
                multi(value, getword(input[i], 4));
                tmp_a = new Capacitor (value, n1 , n2, "C");
                component.push_back(tmp_a);
            } 
            if(getword(input[i], 1)[0] == 'L'){
                multi(value, getword(input[i], 4));
                tmp_a = new Inductor (value, n1, n2, "L");
                component.push_back(tmp_a);
            }
            if(getword(input[i], 1)[0] == 'D'){
                tmp_a = new Diode(n1, n2, "D");
                component.push_back(tmp_a);
            }
            if(getword(input[i], 1)[0] == 'Q'){
                std::string node3;
                int n3 = 0;
                node3 = getword(input[i], 4);
                if(node3 == "0"){
                    n3 = 0;
                }
                else{
                    n3 = stoi(node3.substr(1));
                }
                if(n3 > max_node){
                    max_node = n3;
                }
                std::string type;
                type = getword(input[i], 5);
                tmp_a = new bjt(n1 , n2 , n3 , type, "Q");
                component.push_back(tmp_a);
            }
            if(getword(input[i], 1)[0] == 'M'){
                std::string node3;
                int n3 = 0;
                node3 = getword(input[i], 4);
                if(node3 == "0"){
                    n3 = 0;
                }
                else{
                    n3 = std::stoi(node3.substr(1));
                }
                if(n3 > max_node){
                    max_node = n3;
                }
                std::string type;
                type = getword(input[i], 5);
                tmp_a = new mosfet(n1 , n2 , n3 , type, "M");
                component.push_back(tmp_a);
            }
            if(getword(input[i], 1)[0] == 'G'){
                std::string node3, node4;
                int n3, n4;
                node3 = getword(input[i], 4);
                if(node3 == "0"){
                    n3 = 0;
                }
                else{
                    n3 = std::stoi(node3.substr(1));
                }
                if(n3 > max_node){
                    max_node = n3;
                }
                node4 = getword(input[i], 5);
                if(node4 == "0"){
                    n4 = 0;
                }
                else{
                    n4 = std::stoi(node4.substr(1));
                }
                if(n4 > max_node){
                    max_node = n4;
                }
                multi(value, getword(input[i], 6));
                tmp_a = new v_currsrc(n1, n2, n3, n4, value, "G");
                component.push_back(tmp_a);
            }
        }
    }
}
/*
int main(){
    vector <string> s;
    vector <general*> g;
    s = ReadInput("input3.txt");
    int maxn = 0;
for(int i = 0; i < s.size(); i++){
     std::cout << s[i] << std::endl;
 }

    setting(s,g,maxn);

    for(int i = 0; i < g.size(); i++){
        cout<<g[i]->conductance(100);
    }
}
*/
// int main(){
// vector <string> s;
//    vector <general*> g;
//    s = ReadInput("finaltest.txt");
//     int maxn = 0;
//     setting(s,g,maxn);
//     cout << "testing for setting" << endl;
//     Resistor* tmp2 = dynamic_cast<Resistor*> (g[0]);
//     cout << "designator: " << tmp2->get_type() << "\t" << "value: " << tmp2->resistance
//     << "\t" << "node1: " << tmp2->get_polarity()[0] << "\t" << "node2: " <<
//     tmp2->get_polarity()[1] << "\t" << endl;
//     Resistor* tmp3 = dynamic_cast<Resistor*> (g[1]);
//     cout << "designator: " << tmp3->get_type() << "\t" << "value: " << tmp3->resistance
//     << "\t" << "node1: " << tmp3->get_polarity()[0] << "\t" << "node2: " <<
//     tmp3->get_polarity()[1] << "\t" << endl;
//     Resistor* tmp4 = dynamic_cast<Resistor*> (g[2]);
//     cout << "designator: " << tmp4->get_type() << "\t" << "value: " << tmp4->resistance
//     << "\t" << "node1: " << tmp4->get_polarity()[0] << "\t" << "node2: " <<
//     tmp4->get_polarity()[1] << "\t" << endl;
//     bjt* tmp8 = dynamic_cast<bjt*> (g[3]);
//     cout << "designator: "<< tmp8->get_type() << "\t" << "base node: " << tmp8->b() << "\t"
//     << "collector node: " << tmp8->c() << "\t" << "emitter node: " << tmp8->e() << "\t" <<
//     "type: " << tmp8->type_b << endl;
//     Resistor* tmp5 = dynamic_cast<Resistor*> (g[4]);
//     cout << "designator: " << tmp5->get_type() << "\t" << "value: " << tmp5->resistance
//     << "\t" << "node1: " << tmp5->get_polarity()[0] << "\t" << "node2: " <<
//     tmp5->get_polarity()[1] << "\t" << endl;
//     voltsrc* tmp0 = dynamic_cast<voltsrc*> (g[5]);
//     cout << "designator: " << tmp0->search_name()  << "\t" << "value: " << tmp0->get_volt() 
//     << "\t" << "+ node: " << tmp0->get_polarity()[0] << "\t" << "- node: " << 
//     tmp0->get_polarity()[1] << endl;
//     voltsrc* tmp1 = dynamic_cast<voltsrc*> (g[6]);
//     cout << "designator: " << tmp1->search_name()  << "\t" << "value: " << tmp1->get_volt() 
//     << "\t" << "+ node: " << tmp1->get_polarity()[0] << "\t" << "- node: " << 
//     tmp1->get_polarity()[1] << endl;
//     Capacitor* tmp7 = dynamic_cast<Capacitor*> (g[7]);
//     cout << "designator: " << tmp7->get_type() << "\t" << "value: " << tmp7->capacitance
//     << "\t" << "node1: " << tmp7->get_polarity()[0] << "\t" << "node2: " <<
//     tmp7->get_polarity()[1] << "\t" << endl;
// }

