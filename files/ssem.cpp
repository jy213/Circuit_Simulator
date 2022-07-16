#include<iostream>
#include<complex>
#include<vector>
#include<cmath>
#include<string>
#include<eigen3/Eigen/Dense>
#include"_DCanalysis.cpp"
std::vector<double> frequency (std::vector<double> ac){
    std::vector<double> tmp;
    double n = log10(ac[2] / ac[1]);
    double fs = ac[1];
    double fn = 0;
    for(int i = 0; i < n; i++){
        for(int g = 0; g <= ac[0]; g++ ){
           fn = fs * pow(10, g/ac[0]);
           tmp.push_back(fn);
        }
        fs = 10 * fs;   
    }
    return tmp;
}
void num_of_acvolt ( std::vector<voltsrc*> input, int &num_de, std::vector<int> &g , std::vector<int> &n){
    for (int i = 0; i < input.size(); i++){
        if((input[i]->dc_ac() == "ac") && (input[i]->get_node()[0] != 0)){
            num_de += 1;
            n.push_back(i);
        }
        if((input[i]->dc_ac() == "ac") && (input[i]->get_node()[0] == 0)){
            g.push_back(input[i]->get_node()[1]);
        }
    }
}
std::vector<general*> shortdcsource(std::vector<general*> input, std::vector<shortcircuit>& sc, std::vector<voltsrc*> v){
    std::vector<general*> tmp;
    std::vector<shortcircuit> test;
    shortcircuit val;
    tmp = input;
    general* tmp_str_add;
    for(int i = 0; i < v.size(); i++){
        if(v[i]->dc_ac() == "dc"){
            val.ref_node = v[i]->get_node()[0];
            val.short_node = v[i]->get_node()[1];
            sc.push_back(val);
        }
    }
    for(int m = 0; m < tmp.size(); m++){
        for(int b = 0; b < sc.size(); b++){
            if(tmp[m]->get_polarity()[0] == sc[b].ref_node){
                if(tmp[m]->get_polarity()[1] == sc[b].short_node){
                    tmp.erase(tmp.begin() + m);
                }
            }
            if(tmp[m]->get_polarity()[0] == sc[b].short_node){
                if(tmp[m]->get_polarity()[1] == sc[b].ref_node) {
                    tmp.erase(tmp.begin() + m);
                }
                else{
                    tmp[m]->change_node(sc[b].ref_node, tmp[m]->get_polarity()[1]);
                }
            }
            if(tmp[m]->get_polarity()[1] == sc[b].ref_node){
                if(tmp[m]->get_polarity()[0] == sc[b].short_node){
                    tmp.erase(tmp.begin() + m);
                }
            }
            if(tmp[m]->get_polarity()[1] == sc[b].short_node){
                if(tmp[m]->get_polarity()[0] == sc[b].ref_node){
                    tmp.erase(tmp.begin() + m);
                }
                else{
                    tmp[m]->change_node(tmp[m]->get_polarity()[0], sc[b].ref_node);
                }
            }
        }
    }
    return tmp;
}
void classify_comp_more(std::vector<general*> input, std::vector<Resistor*> &r, std::vector<bjt*> &b, std::vector<mosfet*> &m, std::vector<Diode*>&d, std::vector<voltsrc*> &v , std::vector<currsrc*> &a,std::vector<v_currsrc*>& x){
    for(int i =0; i < input.size(); i++){
        if(input[i]->get_type() == "R"){
            Resistor* tmp = dynamic_cast<Resistor*> (input[i]);
            r.push_back(tmp);
        }
        if(input[i]->get_type() == "Q"){
            bjt* t = dynamic_cast<bjt*> (input[i]);
            b.push_back(t);
        }
        if(input[i]->get_type() == "M"){
            mosfet* t = dynamic_cast<mosfet*> (input[i]);
            m.push_back(t);
        }
        if(input[i]->get_type() == "D"){
            Diode* t = dynamic_cast<Diode*> (input[i]);
            d.push_back(t);
        }
        if(input[i]->get_type() == "V"){

            voltsrc* a = dynamic_cast<voltsrc*> (input[i]);
            v.push_back(a);
        }
        if(input[i]->get_type() == "I"){
            currsrc* m = dynamic_cast<currsrc*> (input[i]);
            a.push_back(m);
        }
        if(input[i]->get_type() == "G"){
            v_currsrc* z = dynamic_cast<v_currsrc*>(input[i]);
            x.push_back(z);
        }
    }
}
Eigen::MatrixXcd SSEM (std::vector<general*> input , int maxnode, double frequency, Eigen::MatrixXd volt, std::vector<shortcircuit> sc){
    std::vector<double> standard_volt;
    for( int e = 0; e < volt.rows(); e++){
        standard_volt.push_back(volt(e,0));
    }
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector<currsrc*> a;
    std::vector<v_currsrc*> tmp_str_voltcurrc; //this vector is used to store the general type pointer to the voltage dependent current source.
    classify_comp_more(input, r, b, m, d, v, a, tmp_str_voltcurrc);
    int num_de = 0;
    std::vector<int> g_volt;
    std::vector<int> n_volt;
    num_of_acvolt(v, num_de, g_volt, n_volt);
    std::vector<int> stand_posi_in_matrix;
    stand_posi_in_matrix = short_str(sc, maxnode);
    int num = num_de + maxnode - sc.size();
    Eigen::MatrixXcd mat(num, num);
    for(int standposi = 0; standposi < num; standposi++){
        int row = 0;
        if (standposi < num_de){
            row = standposi;
        }
        else{
            row = stand_posi_in_matrix[standposi - num_de] + num_de -1;
        }
        std::complex<double> tot_cond (0,0);
        std::vector<int > tmp2;
        tmp2.push_back(0);
        tmp2.push_back(row + 1);
        bool decide = false;       
        for(int c=0; c<num; c++){
            mat(standposi, c) = 0;
        }
        for(int k = 0; k < g_volt.size(); k++){
            if(row == num_de + g_volt[k]-1){
                mat(standposi, find_index(stand_posi_in_matrix , g_volt[k])) = 1;
                decide = true;
            }
        }
        if(row < num_de ){
            int posi = v[n_volt[row]]->get_polarity()[0];
            int nega = v[n_volt[row]]->get_polarity()[1];
            mat(standposi, find_index(stand_posi_in_matrix , posi)) = 1;
            mat(standposi, find_index(stand_posi_in_matrix , nega)) = -1;
            decide = true;                
        }
        std::complex<double> ground_cond (0,0);
        if(!decide){
            for(int colposi = 0; colposi < num; colposi++){
                int index = stand_posi_in_matrix[colposi] -1;
                std::vector<int> tmp;
                if(row - num_de < index){
                    tmp.push_back((row-num_de+1));
                    tmp.push_back((index+1));
                }
                else{
                    tmp.push_back((index+1));
                    tmp.push_back((row-num_de+1));
                }
                std::complex<double> cond (0,0);
                ground_cond = (0,0);
                for(int seq_comp = 0; seq_comp < input.size() ; seq_comp ++){
                    if((input[seq_comp]-> nodes() == 2)&&(input[seq_comp]->get_node()== tmp)){
                        cond = cond - input[seq_comp]-> conductance(frequency * 2 * PI);
                    }
                    if((input[seq_comp]-> nodes() == 2) && (input[seq_comp]->get_node()==tmp2)){
                        ground_cond = ground_cond + input[seq_comp]-> conductance(frequency * 2 * PI);
                    }//this function is used to find volt dependent current source, translate the current output to the conductance at the two control node to the connect node.  
                }
                if(row - num_de != index){
                    mat(standposi, colposi) = cond;
                }
                tot_cond = tot_cond - cond;
            }
            tot_cond += ground_cond;
            mat(standposi, standposi-num_de) = tot_cond;
            if(num_de > 0){
                for (int t = 0; t < num_de ; t++){
                    if(v[n_volt[t]]->get_polarity()[0] == (row - num_de + 1)){
                        mat(standposi, t + maxnode - sc.size()) = -1;
                    }
                    else if(v[n_volt[t]]->get_polarity()[1] == (row - num_de + 1)){
                        mat(standposi , t + maxnode - sc.size()) = 1;
                    }
                }
            }
        }
    }
    for(int h = 0; h < tmp_str_voltcurrc.size(); h++){ //here we add the volt dependent current source
        int v1, v2, c1, c2;
        v1 = find_index(stand_posi_in_matrix, tmp_str_voltcurrc[h]->get_polarity()[0]);
        v2 = find_index(stand_posi_in_matrix, tmp_str_voltcurrc[h]->get_polarity()[1]);
        c1 = find_index(stand_posi_in_matrix, tmp_str_voltcurrc[h]->get_polarity()[2]);
        c2 = find_index(stand_posi_in_matrix, tmp_str_voltcurrc[h]->get_polarity()[3]);
        std::complex<double> trans = tmp_str_voltcurrc[h]->conductance(0);
        if(tmp_str_voltcurrc[h]->get_polarity()[0] != 0){
            if(tmp_str_voltcurrc[h]->get_polarity()[2] != 0){
                mat(v1 + num_de, c1) -= trans;
            }
            if(tmp_str_voltcurrc[h]->get_polarity()[3] != 0){
                mat(v1 +num_de, c2) += trans;
            }
        }
        if(tmp_str_voltcurrc[h]->get_polarity()[1] != 0){
            if(tmp_str_voltcurrc[h]->get_polarity()[2] != 0){
                mat(v2 + num_de, c1) += trans;
            }
            if(tmp_str_voltcurrc[h]->get_polarity()[3] != 0){
                mat(v2 + num_de, c2) -= trans;
            }
        }
    }
    return mat;
}
std::vector<general*> reorganize(std::vector<general*> input , int maxnode, std::vector<double> standard_volt){
    std::vector<general*> tmp_str = input;
    int num_de = standard_volt.size() - maxnode;
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector<currsrc*>a;
    general* tmp;
    classify_comp(input, r, b, m, d, v, a);
    for( int t = 0; t < b.size() ; t ++){
        double VB, VC, VE;
        if(b[t]->e() != 0){
            VE = standard_volt[b[t]->e()-1];
        }
        else{
            VE = 0;
        }
        if(b[t]->b() != 0){
            VB = standard_volt[b[t]->b()-1];
        }
        else{
            VB = 0;
        }
        if(b[t]->c() != 0){
            VC = standard_volt[b[t]->c()-1];
        }
        else{
            VC = 0;
        }
        b[t]->initialize_bjt(VB, VE, VC);
        // store gbe
        tmp = new Resistor(b[t]->get_polarity()[1], b[t]->get_polarity()[2], b[t]->rbe().real(), "R" );
        tmp_str.push_back(tmp);
        // store g0
        tmp = new Resistor(b[t]->get_polarity()[0], b[t]->get_polarity()[2], b[t]->r0().real(), "R");
        tmp_str.push_back(tmp);
        // store gm
        tmp = new v_currsrc(b[t]->get_polarity()[2], b[t]->get_polarity()[0], b[t]->get_polarity()[1], b[t]->get_polarity()[2], b[t]->gm().real() , "G" );
        tmp_str.push_back(tmp);
    }
    for( int p = 0; p < m.size(); p++){
        double VG, VD, VS;
        if(m[p]->d() != 0){
            VD = standard_volt[m[p]->d()-1];
        }
        else{
            VD = 0;
        }
        if(m[p]->s() != 0){
            VS = standard_volt[m[p]->s()-1];
        }
        else{
            VS = 0;
        }
        if(m[p]->g() != 0){
            VG = standard_volt[m[p]->g()-1];
        }
        else{
            VG = 0;
        }
        m[p]->initialize_mos(VG, VS, VD);
        tmp = new Resistor(m[p]->d(), m[p]->s(), m[p]->r0(), "R");
        tmp_str.push_back(tmp);
        tmp = new v_currsrc(m[p]->s(), m[p]->d() , m[p]->g(), m[p]->s(), m[p]->gm(), "G");
        tmp_str.push_back(tmp);
    }
    for(int n = 0; n < d.size(); n++){
        tmp = new Resistor(d[n]->get_polarity()[0], d[n]->get_polarity()[1], d[n]->g(standard_volt[d[n]->get_polarity()[0] + num_de -1], standard_volt[d[n]->get_polarity()[1]+ num_de -1]), "R");
        tmp_str.push_back(tmp);
    }
    return tmp_str;
}
Eigen::MatrixXcd build_acb (std::vector<general*> input, int maxnode, std::vector<shortcircuit> sc){
    int num_de = 0;
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector<currsrc*> a;
    classify_comp(input, r, b, m, d, v, a);
    std::vector<int> g_volt;
    std::vector<int> n_volt;
    num_of_acvolt(v, num_de, g_volt, n_volt);
    std::vector<int> stand_posi_in_matrix;
    stand_posi_in_matrix = short_str(sc, maxnode);
    int num = num_de + maxnode-sc.size();
    Eigen::MatrixXcd col_b (num, 1);
    for(int r = 0; r < num; r++){
        col_b (r, 0) = 0;
    }
    for(int r_n_volt = 0; r_n_volt < n_volt.size(); r_n_volt++){
        col_b (r_n_volt, 0) = v[n_volt[r_n_volt]]->get_volt().real();
    }
    for(int i = 0; i < a.size(); i++){ //to decide whether the current source is connected to the two nodes or one node connecting to the ground.
        if(a[i]->get_node()[0] == 0 && a[i]->dc_ac() == "ac"){
            if(a[i]->get_polarity()[0] != 0){
                col_b (find_index(stand_posi_in_matrix ,a[i]->get_polarity()[0]) + num_de, 0) = a[i]->get_current().real();
            }
            else{
                col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[1])+ num_de, 0) = - a[i]->get_current().real();
            }
        }
        if(a[i]->get_node()[0] != 0 && a[i]->dc_ac() == "ac"){
            col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[0])+ num_de, 0) = a[i]->get_current().real();
            col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[1])+ num_de, 0) = - a[i]->get_current().real();
        }
    }
    for(int i = 0; i < v.size(); i++){
        if(v[i]->get_node()[0] == 0){
            if(v[i]->get_polarity()[0] != 0){
                col_b (find_index(stand_posi_in_matrix, v[i]->get_node()[1]) + num_de, 0) = v[i]->get_volt().real();
            }
            if(v[i]->get_polarity()[1] != 0){
                col_b (find_index(stand_posi_in_matrix, v[i]->get_node()[1]) + num_de, 0) = - v[i]->get_volt().real();
            }
        }
    }
    return col_b;
}
Eigen::MatrixXcd recover_complex_circuit(Eigen::MatrixXcd stand_volt, std::vector<shortcircuit> sc){
    Eigen::MatrixXcd tmp (stand_volt.rows() + sc.size(), 1);
    for(int t = 0; t < tmp.rows(); t++){
        tmp(t, 0) = 0;
    }
    std::vector<int> node_in_short = short_str(sc, tmp.rows());
    for(int g = 0; g < node_in_short.size(); g++){ //fill in gap that do not need alter position
        tmp(node_in_short[g] -1 , 0) = stand_volt(g,0);
    }
    for(int m = 0; m < sc.size(); m++){ //fill the shorted one
        if(sc[m].ref_node == 0){
            tmp(sc[m].short_node -1, 0) = 0;
        }
        else{
        tmp(sc[m].short_node -1, 0) = tmp(sc[m].ref_node-1);
        }
    }
    return tmp;
}
void find_final_sol (){
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector<currsrc*> a;
    std::cout<< "please enter the name of text file as input: " <<std::endl;
    std::string file;
    std::cin >> file;
    std::vector <std::string> s;
    s = ReadInput(file);
    std::vector <general*> in;
    int maxnode = 0;
    setting(s, in, maxnode);
    std::vector<double> f;
    f = frequency(ac(s));
    Eigen::MatrixXd dcvolt;
    dcvolt = get_standart_volt(in, maxnode);
    std::vector<double> tmp;
    for(int i = 0; i < dcvolt.rows(); i++){
        tmp.push_back(dcvolt(i,0));
    }
    std::complex<double> inp = (0,0);
    std::vector <std::string> insrc;
    std::vector<voltsrc*> vstore;
    std::vector<currsrc*> istore;
    for(int j=0; j<in.size(); j++){
        if(in[j]->get_type() == "V"){
            voltsrc* a = dynamic_cast<voltsrc*> (in[j]);
            if(a->dc_ac() == "ac"){
                inp = a->get_volt();
                std::string aa = a->search_name();
                insrc.push_back(aa);
                vstore.push_back(a);
            }
        }
        if(in[j]->get_type() == "I"){
            currsrc* b = dynamic_cast<currsrc*> (in[j]);
            if(b->dc_ac() == "ac"){
                inp = b->get_current();
                std::string bb = b->search_name();
                insrc.push_back(bb);
                istore.push_back(b);
            }
        }      
    }
    int ins = -1;
    if(insrc.size() > 1){
        std::cout<< "multiple sources detected, please choose one from the lists below as input: " <<std::endl;
        for(int q=0; q<vstore.size(); q++){
            std::cout<< q << ": " << vstore[q]->search_name() << std::endl;
        }
        for(int p=vstore.size(); p<vstore.size()+istore.size(); p++){
            std::cout<< p << ": " << istore[p]->search_name() << std::endl;
        }
        std::cin >> ins; 
        if(0  <= ins < vstore.size()){
            // volt source
            inp = vstore[ins]->get_volt();
        }
        if( vstore.size()<=ins) {
            // current source
            inp = istore[ins - vstore.size()]->get_current();
        }
    }
    std::cout<< "please choose a node as output: (type a number) "<< std::endl;
    std::cout<< "choose between 1 ~ " << maxnode << std::endl;
    int node;
    std::cin >> node;
    std::vector<general*> testn;
    setting(s, testn, maxnode);
    std::vector<general*> input = reorganize(testn , maxnode, tmp);
    std::vector <double> mag;
    std::vector <double> phase;
    classify_comp(input, r, b, m, d, v, a);
    std::vector<shortcircuit> shortc;
    std::vector<general*> finalinput = shortdcsource(input ,shortc ,v);
    for(int y = 0; y < f.size(); y++){
        Eigen::MatrixXcd G = SSEM(finalinput, maxnode, f[y], dcvolt, shortc);
        Eigen::MatrixXcd colb = build_acb(finalinput, maxnode, shortc);
        Eigen::MatrixXcd result = G.colPivHouseholderQr().solve(colb);
        Eigen::MatrixXcd final = recover_complex_circuit(result, shortc);
        std::complex<double> trans;
        trans.real(final(node-1,0).real());
        trans.imag(final(node-1,0).imag());
        //calculate magnitude and phase of transfer function
        mag.push_back(sqrt(pow(trans.real(), 2) + pow(trans.imag(), 2))/sqrt(pow(inp.real(), 2) + pow(inp.imag() , 2)));
        if(trans.real()<0 && trans.imag()<0){
            phase.push_back(180 * atan(trans.imag() / trans.real()) / PI -180);
        }
        if(trans.real()<0 && trans.imag()>0){
            phase.push_back(180 * atan(trans.imag() / trans.real()) / PI + 180);
        }
        if(trans.real()>0){
            phase.push_back(180 * atan(trans.imag() / trans.real()) / PI);
        }
        if(trans.real()==0){
            phase.push_back(0);
        }
    }
    std::ofstream outFile;
    std::cout<< "please enter the name of the text file to store the output: " << std::endl;
    std::string ot;
    std::cin >> ot;
    outFile.open(ot);
    outFile << "frequency" << "\t" << "magnitude" << "\t" << "phase" << std::endl;
    for(int k=0;k<f.size();k++){
      outFile << f[k] << "\t" << mag[k] << "\t" << phase[k] << std::endl;  
    }
    outFile.close();
    std::cout<< "completed!" << std::endl;
}
//int main(){
  //  find_final_sol();
//}
/*
int main(){       
    vector <string> s;
    vector <general*> g;
    s = ReadInput("mosfet.txt");
    int maxn = 0;
    setting(s,g,maxn);
    vector<general*> testn;
    Eigen::MatrixXd dcvolt;
    dcvolt = get_standart_volt(g, maxn);
    std::vector<double> tmp;
    for( int i = 0; i < dcvolt.rows(); i++){
        tmp.push_back(dcvolt(i,0));
    }
    setting(s, testn, maxn);
    std::vector<general*> test = reorganize(testn, maxn, tmp);
    int num_de = 0;
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector<currsrc*> a;
    classify_comp(test, r, b, m, d, v, a);
    std::vector<shortcircuit> test3;
    std::vector<general*> test2 = shortdcsource(test,test3 ,v);
    Eigen::MatrixXcd testcol = SSEM(test2, maxn, 100, dcvolt, test3);
    Eigen::MatrixXcd test5 = build_acb(test2, maxn, test3);
    std::cout << "result for build_acb: " << std::endl;
    std::cout << test5 << std::endl;
    std::cout << std::endl;
    std::cout << "result for SSEM: " << std::endl;
    std::cout << testcol<< std::endl;
    std::cout << std::endl;
    Eigen::MatrixXcd result = testcol.colPivHouseholderQr().solve(test5);
    Eigen::MatrixXcd final = recover_complex_circuit(result, test3);
    std::cout << "result for recover_complex_circuit: " << std::endl;
    std::cout << final << std::endl;
}
*/