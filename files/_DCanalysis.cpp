#include<iostream>
#include<complex>
#include<vector>
#include<cmath>
#include<string>
#include"part1.cpp"
#include<eigen3/Eigen/Dense>
struct shortcircuit{
    int ref_node;
    int short_node;
};
std::vector<general*> short_circuit(std::vector<general*> in, std::vector<shortcircuit> &str_posi){
    std::vector<general*> tmp;
    std::vector<general*> new_tmp;
    std::vector<shortcircuit> test;
    shortcircuit val;
    tmp = in;
    general* tmp_str_add;
    for(int i = 0; i < tmp.size(); i++){
        if(tmp[i]->get_type() == "L"){
            val.ref_node = tmp[i]->get_node()[0];
            val.short_node = tmp[i]->get_node()[1];
            test.push_back(val);
            str_posi.push_back(val);
        }
        if(tmp[i]->get_type() == "V"){
            voltsrc* tmpc = dynamic_cast<voltsrc*> (tmp[i]);
            if(tmpc->dc_ac() == "ac") {
                val.ref_node = tmp[i]->get_node()[0];
                val.short_node = tmp[i]->get_node()[1];
                test.push_back(val);
                str_posi.push_back(val);
            }
        }
    }
    for(int m = 0; m < tmp.size(); m++){
        if(tmp[m]->nodes() == 2 || tmp[m]->nodes() == 4){
            for(int b = 0; b < test.size(); b++){
                if(tmp[m]->get_polarity()[0] == test[b].ref_node){
                    if(tmp[m]->get_polarity()[1] == test[b].short_node){
                        tmp.erase(tmp.begin() + m);
                    }
                }
                if(tmp[m]->get_polarity()[0] == test[b].short_node){
                    if(tmp[m]->get_polarity()[1] == test[b].ref_node) {
                        tmp.erase(tmp.begin() + m);
                    }
                    else{
                        tmp[m]->change_node(test[b].ref_node, tmp[m]->get_polarity()[1]);
                    }
                }
                if(tmp[m]->get_polarity()[1] == test[b].ref_node){
                    if(tmp[m]->get_polarity()[0] == test[b].short_node){
                        tmp.erase(tmp.begin() + m);
                    }
                }
                if(tmp[m]->get_polarity()[1] == test[b].short_node){
                    if(tmp[m]->get_polarity()[0] == test[b].ref_node) {
                        tmp.erase(tmp.begin() + m);
                    }
                    else{
                        tmp[m]->change_node(tmp[m]->get_polarity()[1], test[b].ref_node);
                    }
                }
            }
        }
        if(tmp[m]->nodes() == 3){
            if(tmp[m]->get_type() == "Q"){
                bjt* tmpb = dynamic_cast<bjt*> (tmp[m]);
                for( int c = 0; c < test.size() ; c++){
                    if(tmpb->get_polarity()[0] == test[c].short_node){
                        tmpb->change_node(test[c].ref_node, tmpb->get_polarity()[1], tmpb->get_polarity()[2]);
                    }
                    if(tmpb->get_polarity()[1] == test[c].short_node){
                        tmpb->change_node(tmpb->get_polarity()[0], test[c].ref_node, tmpb->get_polarity()[2]);
                    }
                    if(tmpb->get_polarity()[2] == test[c].short_node){
                        tmpb->change_node(tmpb->get_polarity()[0], tmpb->get_polarity()[1],test[c].ref_node);
                    }
                }
                tmp.erase(tmp.begin() + m);
                tmp_str_add = new bjt(tmpb->get_polarity()[0], tmpb->get_polarity()[1], tmpb->get_polarity()[2], tmpb->getmodel(), tmpb->get_type());
                tmp.push_back(tmp_str_add);
            }
            if(tmp[m]->get_type() == "M"){
                mosfet* tmpb = dynamic_cast<mosfet*> (tmp[m]);
                    for( int c = 0; c < test.size() ; c++){
                    if(tmpb->get_polarity()[0] == test[c].short_node){
                        tmpb->change_node(test[c].ref_node, tmpb->get_polarity()[1], tmpb->get_polarity()[2]);
                    }
                    if(tmpb->get_polarity()[1] == test[c].short_node){
                        tmpb->change_node(tmpb->get_polarity()[0],test[c].ref_node, tmpb->get_polarity()[2]);
                    }
                    if(tmpb->get_polarity()[2] == test[c].short_node){
                        tmpb->change_node(tmpb->get_polarity()[0], tmpb->get_polarity()[1],test[c].ref_node);
                    }
                }
                tmp.erase(tmp.begin() + m);
                tmp_str_add = new mosfet(tmpb->d(), tmpb->g(), tmpb->s(), tmpb->getmodel(), tmpb->get_type());
                tmp.push_back(tmp_str_add);
            }
        }
    }
    return tmp;
}
std::vector<int> short_str(std::vector<shortcircuit> sc, int maxnodein){//this function is used to classify the ref nodes and put them into a vector
    std::vector<int> node_in_short;
    for( int y = 1; y < maxnodein + 1; y++){
        bool decide = true;
        for(int x = 0; x < sc.size(); x++){
            if(y == sc[x].short_node){
                decide = false;
            }
        }
        if(decide){
            node_in_short.push_back(y);
        }  
    }
    return node_in_short;
}
int find_index(std::vector<int> input, int val){
    for(int i = 0; i< input.size(); i++){
        if(input[i] == val){
            return i;
        }
    }
    return 0;
}

void dc_volt (std::vector<voltsrc*> tmp, int &num, std::vector<int> &n, std::vector<int> &g, std::vector<int> &ac, std::vector<int>& vg){
    for (int i = 0; i < tmp.size(); i++){
        if(tmp[i]->get_node()[0] != 0 && tmp[i]->get_node()[1] != 0 && tmp[i]->dc_ac() == "dc"){
            num += 1;
            n.push_back(i);
        }
        if(tmp[i]->get_node()[0] == 0 && tmp[i]->dc_ac() == "dc"){
            g.push_back(tmp[i]->get_node()[1]);
            vg.push_back(i);
        }
        if(tmp[i]->dc_ac() == "ac"){
            ac.push_back(i);
        }
    }        
}
void classify_comp(std::vector<general*> input, std::vector<Resistor*> &r, std::vector<bjt*> &b, std::vector<mosfet*> &m, std::vector<Diode*>&d, std::vector<voltsrc*> &v , std::vector<currsrc*> &a){
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
    }
}
std::vector<general*> reorganizedc(std::vector<general*> input, int& extra_node){
    std::vector<general*> tmp = input;
    voltsrc* tmpv;
    for(int i = 0; i < input.size(); i++){
        if(input[i]->get_type() == "Q"){
            bjt* tmp_bjt = dynamic_cast<bjt*> (input[i]);
            if(!(tmp_bjt->b()==0 && tmp_bjt->e()== 0)){
            tmpv = new voltsrc(0.7, tmp_bjt->get_polarity()[1] , tmp_bjt->get_polarity()[2], "V" ,"vn");
            tmp.push_back(tmpv);
            }
            if(tmp_bjt->get_polarity()[1]!= 0 && tmp_bjt->get_polarity()[2]!= 0 ){
                extra_node +=1;
            }
            if(!(tmp_bjt->c()==0 && tmp_bjt->e() == 0)){
            tmpv = new voltsrc(1, tmp_bjt->get_polarity()[0] , tmp_bjt->get_polarity()[2], "V", "vn");
            tmp.push_back(tmpv);
            }
            if(tmp_bjt->get_polarity()[0]!= 0 && tmp_bjt->get_polarity()[2]!= 0 ){
                extra_node +=1;
            }
        }
        if(input[i]->get_type() == "M"){
            mosfet * tmpm = dynamic_cast<mosfet*> (input[i]);
            if(!(tmpm->g()==0 &&tmpm->s()==0)){
            voltsrc* tmpv = new voltsrc(1.1 * abs(tmpm->vt()), tmpm->g() , tmpm->s(), "V", "vn"); 
            tmp.push_back(tmpv);
            }
            if(tmpm->g()!=0 && tmpm->s()!= 0 ){
                extra_node += 1;
            }
            if(!(tmpm->d() == 0 && tmpm->s() == 0)){
            tmpv = new voltsrc(2* abs(tmpm->vt()), tmpm->d() , tmpm->s(), "V", "vn");
            tmp.push_back(tmpv); 
            }
        }
        if(input[i]->get_type() == "D"){
            Diode* tmpd = dynamic_cast<Diode*>(input[i]);
            voltsrc* tmpv = new voltsrc(0.7, tmpd->get_polarity()[0] , tmpd->get_polarity()[1], "V", "vn"); 
            tmp.push_back(tmpv);
            if(tmpd->get_polarity()[0]!= 0 && tmpd->get_polarity()[1]!= 0){
                    extra_node +=1;
            }   
        }
    }
    return tmp;
}
Eigen::MatrixXd recover_circuit(Eigen::MatrixXd stand_volt, std::vector<shortcircuit> sc){
    Eigen::MatrixXd tmp (stand_volt.rows() + sc.size(), 1);
    for(int t = 0; t < tmp.rows(); t++){
        tmp(t, 0) = 0;
    }
    std::vector<int> node_in_short = short_str(sc, tmp.rows());
    // fill in gap that do not need alter position
    for(int g = 0; g < node_in_short.size(); g++){
        tmp(node_in_short[g] -1 , 0) = stand_volt(g,0);
    }
    // fill the shorted ones
    for(int m = 0; m < sc.size(); m++){
        if(sc[m].ref_node == 0){
            tmp(sc[m].short_node -1, 0) = 0;
        }
        else{
        tmp(sc[m].short_node -1, 0) = tmp(sc[m].ref_node-1);
        }
    }
    return tmp;
}
//the three functions below are the core of Newton Raphson matrix
Eigen::MatrixXd build_iterate_matrix (std::vector<general*> in, int maxnode, Eigen::MatrixXd guess){
    std::vector<shortcircuit> shortc;
    std::vector<general*> input = short_circuit(in, shortc);
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector<double> volt_vector;
    std::vector<currsrc*> a;
    std::vector<v_currsrc*> vc; 
    for(int i = 0; i < guess.rows(); i++){
        volt_vector.push_back(guess(i,0));
    }
    classify_comp(input, r, b, m, d, v, a);
    int num_de = 0;
    std::vector<int> g_volt;
    std::vector<int> n_volt; //this vector is used to store the nodes of voltsrc connecting to ground.
    std::vector<int> v_ground_posi; //this vector contains the index for the voltsrc connecting to two nodes.
    std::vector<int> v_ac;
    std::vector<int> stand_posi_in_matrix;
    stand_posi_in_matrix = short_str(shortc, maxnode);
    dc_volt(v, num_de, n_volt, g_volt,v_ac, v_ground_posi);
    int num = num_de + maxnode - shortc.size();
    Eigen::MatrixXd mat(num, num);
    for(int standposi = 0; standposi < num; standposi++){
        int row = 0;
        if (standposi < num_de){
            row = standposi;
        }
        else{
            row = stand_posi_in_matrix[standposi - num_de] + num_de -1;
        }
        double tot_cond= 0;
        bool decide = false;
        for(int c=0; c<num; c++){
            mat(standposi, c) = 0;
        }
        // construct model row full 0;
        for(int k = 0; k < g_volt.size(); k++){
            if(row == num_de + g_volt[k]-1){ //connect to ground, put row, attention to the position of this.
                mat(standposi, find_index(stand_posi_in_matrix , g_volt[k])) = 1;
                decide = true;
            }
        }
        if(row < num_de ){ // connect to two nodes
            int posi = v[n_volt[row]]->get_polarity()[0];
            int nega = v[n_volt[row]]->get_polarity()[1];
            mat(standposi, find_index(stand_posi_in_matrix , posi)) = 1;
            mat(standposi, find_index(stand_posi_in_matrix , nega)) = -1;
            decide = true;                
        }
        double ground_cond= 0;
        //construction of two possible voltages completed, now -> components
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
                std::vector<int > tmp2; //construct compare vector;
                tmp2.push_back(0);
                tmp2.push_back(row + 1);
                double cond = 0;
                ground_cond =0;
                vc.clear();
                for(int seq_comp = 0; seq_comp < input.size() ; seq_comp ++){
                    if((input[seq_comp]->get_node()== tmp)&&(input[seq_comp]->get_type() == "R")){
                        cond = cond - input[seq_comp]-> conductance(0).real();
                    }
                    if((input[seq_comp]->get_node()==tmp2) && (input[seq_comp]->get_type() == "R")){
                        ground_cond = ground_cond + input[seq_comp]-> conductance(0).real();
                    }
                    //this function is used to find volt dependent current source.
                    //translate the current output to the conductance at the two control node to the connect node.
                    if(input[seq_comp]->get_type() == "G") { 
                        v_currsrc* mid = dynamic_cast<v_currsrc*> (input[seq_comp]);
                        vc.push_back(mid); 
                    }
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
                        mat(standposi, t + maxnode - shortc.size()) = -1;
                    }
                    else if(v[n_volt[t]]->get_polarity()[1] == (row - num_de + 1)){
                        mat(standposi , t + maxnode - shortc.size()) = 1;
                    }
                }
            }
        }
    }
    for(int h = 0; h < vc.size(); h++){
        //here we add the volt dependent current source
        int v1, v2, c1, c2;
        v1 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[0]);
        v2 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[1]);
        c1 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[2]);
        c2 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[3]);
        double trans = vc[h]->conductance(0).real();
        if(vc[h]->get_polarity()[0] != 0){
            if(vc[h]->get_polarity()[2] != 0){
                mat(v1 + num_de, c1) -= trans;
            }
            if(vc[h]->get_polarity()[3] != 0){
                mat(v1 +num_de, c2) += trans;
            }
        }
        if(vc[h]->get_polarity()[1] != 0){
            if(vc[h]->get_polarity()[2] != 0){
                mat(v2 + num_de, c1) += trans;
            }
            if(vc[h]->get_polarity()[3] != 0){
                mat(v2 + num_de, c2) -= trans;
            }
        }
    }
    for (int i = 0; i < d.size(); i++){ //for diode      
        int anode = find_index(stand_posi_in_matrix , d[i]->get_polarity()[0]);
        int cathode = find_index(stand_posi_in_matrix , d[i]->get_polarity()[1]);
        if(d[i]->get_node()[0]!= 0){
            //anode row
            mat(anode + num_de, cathode) += - d[i]->g(volt_vector[anode], volt_vector[cathode]);
            mat(anode + num_de, anode) += d[i]->g(volt_vector[anode], volt_vector[cathode]);
            //cathode row
            mat(cathode + num_de, cathode) += d[i]->g(volt_vector[anode], volt_vector[cathode]);
            mat(cathode + num_de, anode) += -d[i]->g(volt_vector[anode], volt_vector[cathode]);
        }
        else{
            if(d[i]->get_polarity()[0]== 0){
                mat(cathode + num_de, cathode) += d[i]->g(0, volt_vector[cathode]);
            }
            if(d[i]->get_polarity()[1] == 0){
                mat(anode + num_de, anode) += d[i]->g(volt_vector[anode], 0);
            }
        }
    }
    for(int m = 0; m < b.size(); m++){ //for bjt
        double VB, VC, VE;
        int b_node, c_node, e_node;
        c_node = find_index(stand_posi_in_matrix, b[m]->get_polarity()[0]);
        b_node = find_index(stand_posi_in_matrix, b[m]->get_polarity()[1]);
        e_node = find_index(stand_posi_in_matrix, b[m]->get_polarity()[2]);
        if(b[m]->get_polarity()[2] != 0){
             VE = volt_vector[e_node];
        }
        else{
            VE = 0;
        }
        if(b[m]->get_polarity()[1] != 0){
            VB = volt_vector[b_node];
        }
        else{
            VB = 0;
        }
        if(b[m]->get_polarity()[0] != 0){
            VC = volt_vector[c_node];
        }
        else{
            VC = 0;
        }
        b[m]->initialize_bjt(VB, VE, VC);
        if(b[m]->b() != 0){ //for bjt with b on that row.
            mat(b_node + num_de, b_node) += b[m]->b_b();
            if(b[m]->e() != 0){ 
                mat(b_node + num_de, e_node) += b[m]->b_e();
            }
            if(b[m]->c() != 0){
                mat(b_node + num_de, c_node -1) += b[m]->b_c();            
            }
        }
        if(b[m]->e()!= 0){ //for bjt with e on that row.
            mat(e_node + num_de, e_node) -= b[m]->e_e(); 
            if(b[m]->b()!= 0){
                mat(e_node + num_de, b_node) -= b[m]->e_b();
            }
            if(b[m]->c()!=0){
                mat(e_node, c_node) -= b[m]->e_c();
            }   
        }
        if(b[m]->c()!= 0){ //for bjt with c on that row
            if(b[m]->e()!= 0){
                mat(c_node, e_node) += b[m]->c_e();
            }
            if(b[m]->b() != 0){ 
                mat(c_node, b_node) += b[m]->c_b();
            }
            mat(c_node, c_node) += b[m]->c_c();
        }
    }
    for(int g = 0; g < m.size(); g++){ //for mosfet
        double VG, VD, VS;
        int g_node, d_node, s_node;
        g_node = find_index(stand_posi_in_matrix, m[g]->g());
        d_node = find_index(stand_posi_in_matrix, m[g]->d());
        s_node = find_index(stand_posi_in_matrix, m[g]->s());
        if(m[g]->d() != 0){
            VD = volt_vector[d_node];
        }
        else{
            VD = 0;
        }
        if(m[g]->s() != 0){
            VS = volt_vector[s_node];
        }
        else{
            VS = 0;
        }
        if(m[g]->g() != 0){
            VG = volt_vector[g_node];
        }
        else{
            VG = 0;
        }
        m[g]->initialize_mos(VG, VS, VD);
        if(m[g]->s() != 0){ //for mos with s on that row
            if(m[g]->g()!= 0){
                mat(s_node + num_de, g_node) -= m[g]->Id_g();
            }
            mat(s_node + num_de, s_node) -= m[g]->Id_s();
            if(m[g]->d() != 0 ){
                mat(s_node + num_de, d_node) -= m[g]->Id_d();
            }
        }
       if (m[g]->d()!= 0){ //for mos with d on that row
            if(m[g]->g()!= 0){
                mat(d_node + num_de, g_node) += m[g]->Id_g();
            }
            if(m[g]->s() != 0){
                mat(d_node + num_de, s_node) += m[g]->Id_s();
            }
            mat(d_node + num_de, d_node) += m[g]->Id_d();
        }
    }
    return mat;
}
Eigen::MatrixXd build_fvm_matrix(std::vector<general*> in, int maxnode, Eigen::MatrixXd vin){
// to build fvm for each node. It should contains the linear and non-linear terms , the sum of it should be equal to zero for the function. 
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector <double> guess_volt;
    std::vector<currsrc*> a;
    std::vector<shortcircuit> shortc;
    std::vector<general*> input = short_circuit(in, shortc);
    std::vector<int> stand_posi_in_matrix;
    stand_posi_in_matrix = short_str(shortc, maxnode);
    for( int u = 0; u < vin.rows(); u++){
        guess_volt.push_back(vin(u , 0));
    }
    classify_comp(input, r, b, m, d, v, a);
    int num_de = 0;
    std::vector<int> g_volt; //this vector is used to store the nodes of voltsrc connecting to ground.
    std::vector<int> n_volt; //this vector contains the index for the voltsrc connecting to two nodes.
    std::vector<int> ground_p;
    std::vector<int> v_ac;
    dc_volt(v, num_de, n_volt, g_volt, v_ac, ground_p);
    int num = num_de + maxnode - shortc.size();
    Eigen::MatrixXd col_b (num, 1);
    for(int r = 0; r < num; r++){
        col_b (r, 0) = 0;
    }
    for(int r_n_volt = 0; r_n_volt < n_volt.size(); r_n_volt++){ //want vposi - vnega -vsrc = 0
        col_b (r_n_volt, 0) = guess_volt[ find_index(stand_posi_in_matrix,v[n_volt[r_n_volt]]->get_polarity()[0])] - guess_volt[find_index(stand_posi_in_matrix,v[n_volt[r_n_volt]]->get_polarity()[1])] - v[n_volt[r_n_volt]]->get_volt().real();
    }
    for(int i = 0; i < a.size(); i++){
        if(a[i]->get_node()[0] == 0 && a[i]->dc_ac() == "dc"){
            if(a[i]->get_polarity()[0] != 0){
                col_b (find_index(stand_posi_in_matrix ,a[i]->get_polarity()[0]) + num_de, 0) -= a[i]->get_current().real();
            }
            else{
                col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[1])+ num_de, 0) += a[i]->get_current().real();
            }
        }
        if(a[i]->get_node()[0] != 0 && a[i]->dc_ac() == "dc"){
            col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[0])+ num_de, 0) -= a[i]->get_current().real();
            col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[1])+ num_de, 0) += a[i]->get_current().real();
        }
    }
    for(int g = 0; g < d.size(); g++){ //consider diode
        int anode = find_index(stand_posi_in_matrix , d[g]->get_polarity()[0]);
        int cathode = find_index(stand_posi_in_matrix , d[g]->get_polarity()[1]);
        if(d[g]->get_node()[0]!=0){
        col_b (anode + num_de, 0) += d[g]->id(guess_volt[anode + num_de], guess_volt[cathode + num_de]);
        col_b (cathode + num_de,0) -= d[g]->id(guess_volt[anode + num_de], guess_volt[cathode + num_de]);
        }
        if(d[g]->get_polarity()[0] == 0){
            col_b (cathode + num_de, 0) += d[g]->id(0, guess_volt[cathode + num_de]);
        }
        if(d[g]->get_polarity()[1] == 0){
            col_b (anode + num_de, 0) += d[g]->id(guess_volt[anode + num_de], 0);
        }
    }
    for(int h = 0; h < b.size(); h++){ //consider bjt
        double VB, VC, VE;
        int b_node, c_node, e_node;
        c_node = find_index(stand_posi_in_matrix, b[h]->get_polarity()[0]);
        b_node = find_index(stand_posi_in_matrix, b[h]->get_polarity()[1]);
        e_node = find_index(stand_posi_in_matrix, b[h]->get_polarity()[2]);
        if(b[h]->get_polarity()[2] != 0){
            VE = guess_volt[e_node];
        }
        else{
            VE = 0;
        }
        if(b[h]->get_polarity()[1] != 0){
            VB = guess_volt[b_node];
        }
        else{
            VB = 0;
        }
        if(b[h]->get_polarity()[0] != 0){
            VC = guess_volt[c_node];
        }
        else{
            VC = 0;
        }
        b[h]->initialize_bjt(VB, VE, VC);
        if(b[h]->c() != 0){
            col_b(c_node + num_de, 0) += b[h]->ic();
        }
        if(b[h]->b() != 0){ //node B
            col_b(b_node+ num_de, 0) += b[h]->ib();
        }
        if(b[h]->e() != 0){ //node E
            col_b(e_node+ num_de, 0) -= b[h]->ie();
        }
    } 
    for(int u =0; u < m.size(); u ++){ //consider mosfet
        double VG, VD, VS;
        int g_node, d_node, s_node;
        g_node = find_index(stand_posi_in_matrix, m[u]->g());
        d_node = find_index(stand_posi_in_matrix, m[u]->d());
        s_node = find_index(stand_posi_in_matrix, m[u]->s());
        if(m[u]->get_polarity()[2] != 0){
            VD = guess_volt[d_node];
        }
        else{
            VD = 0;
        }
        if(m[u]->get_polarity()[1] != 0){
            VS = guess_volt[s_node];
        }
        else{
            VS = 0;
        }
        if(m[u]->get_polarity()[0] != 0){
            VG = guess_volt[g_node];
        }
        else{
            VG = 0;
        }
        m[u]->initialize_mos(VG, VS, VD);
        if(m[u]->d() != 0){ //node D
            col_b(d_node + num_de, 0) += m[u]->id();
        }
        if(m[u]->s() != 0){ //node S 
            col_b(s_node + num_de, 0) -= m[u]->id();
        }
    }
    for(int y = 0; y < r.size(); y++){ //consider resistor
        int posi, nega;
        posi = find_index(stand_posi_in_matrix, r[y]->get_polarity()[0]);
        nega = find_index(stand_posi_in_matrix, r[y]->get_polarity()[1]);
        if(r[y]->get_node()[0] != 0){
            col_b(posi+ num_de , 0) += (guess_volt[posi] - guess_volt[nega]) * r[y]->conductance(0).real();
            col_b(nega+ num_de , 0) += (guess_volt[nega] - guess_volt[posi]) * r[y]->conductance(0).real();
        }
        else{
            col_b(find_index(stand_posi_in_matrix, r[y]->get_node()[1]) + num_de, 0) += guess_volt[find_index(stand_posi_in_matrix, r[y]->get_node()[1])] * r[y]->conductance(0).real();
        }
    }
    for( int f = 0 ; f < ground_p.size(); f++ ){ //considering the voltage connecting to the ground
        int posiv, negav;
        posiv = find_index(stand_posi_in_matrix, v[ground_p[f]]->get_polarity()[0]);
        negav = find_index(stand_posi_in_matrix, v[ground_p[f]]->get_polarity()[1]);
        if((v[f]->get_type() == "V") && (v[f]->get_node()[0] == 0)){
            col_b (posiv + num_de, 0) =  guess_volt[posiv]-v[ground_p[f]]->get_volt().real();
            }
        if(v[f]->get_polarity()[1] != 0){
            col_b (negav + num_de, 0) = guess_volt[negav] - v[ground_p[f]]->get_volt().real();
        }
    }
    return col_b;
}
Eigen::MatrixXd build_guess_volt (std::vector<general*> in, int maxnode){
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    std::vector<currsrc*> a;
    int extra_node = 0;
    std::vector<general*> input;
    std::vector<shortcircuit> sc;
    std::vector<v_currsrc*> vc;
    input = short_circuit(in, sc);
    input = reorganizedc(input, extra_node);
    classify_comp(input, r, b, m, d, v, a );
    int num_de = 0;
    std::vector<int> g_volt;
    std::vector<int> n_volt;
    std::vector<int> v_ac;
    std::vector<int> ground_p;
    std::vector<int> stand_posi_in_matrix;
    stand_posi_in_matrix = short_str(sc, maxnode);
    dc_volt(v, num_de, n_volt, g_volt, v_ac, ground_p);
    int num = num_de + maxnode - sc.size();
    Eigen::MatrixXd mat(num, num);
    for(int stand_posi = 0; stand_posi < num; stand_posi++){
        int row = 0;
        if (stand_posi < num_de){
            row = stand_posi;
        }
        else{
            row = stand_posi_in_matrix[stand_posi - num_de] + num_de -1;
        }
        std::complex<double> tot_cond (0,0);
        std::vector<int > tmp2;
        tmp2.push_back(0);
        tmp2.push_back(row + 1);
        bool decide = false;       
        for(int c=0; c<num; c++){
            mat(stand_posi, c) = 0;
        }
        for(int k = 0; k < g_volt.size(); k++){
            if(row == num_de + g_volt[k]-1){
                mat(stand_posi, find_index(stand_posi_in_matrix, g_volt[k])) = 1;
                decide = true;
            }
        }
        if(row < num_de ){
            int posi = v[n_volt[row]]->get_polarity()[0];
            int nega = v[n_volt[row]]->get_polarity()[1];
            mat(stand_posi, find_index(stand_posi_in_matrix , posi)) = 1;
            mat(stand_posi, find_index(stand_posi_in_matrix , nega)) = -1;
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
                    if((input[seq_comp]->get_node()== tmp)&&(input[seq_comp]->get_type() == "R")){
                        cond = cond - input[seq_comp]-> conductance(0);
                    }
                    if((input[seq_comp]->get_node()==tmp2) && (input[seq_comp]->get_type() == "R")){
                        ground_cond = ground_cond + input[seq_comp]-> conductance(0);
                    }
                    //this function is used to find volt dependent current source.
                    //translate the current output to the conductance at the two control node to the connect node.  
                    if(input[seq_comp]->get_type() == "G") { 
                        v_currsrc* mid = dynamic_cast<v_currsrc*> (input[seq_comp]);
                        vc.push_back(mid);                           
                    }
                    if(row - num_de != index){
                        mat(stand_posi, colposi) = cond.real();
                    }
                    tot_cond = tot_cond - cond;
                }
                tot_cond += ground_cond;
                mat(stand_posi, stand_posi-num_de) = tot_cond.real();
                if(num_de > 0){
                    for (int t = 0; t < num_de ; t++){
                        if(v[n_volt[t]]->get_polarity()[0] == (row - num_de + 1)){
                            mat(stand_posi, t + maxnode - sc.size()) = -1;
                        }
                        else if(v[n_volt[t]]->get_polarity()[1] == (row - num_de + 1)){
                            mat(stand_posi , t + maxnode - sc.size()) = 1;
                        }
                    }
                }
            }
        }
        for(int h = 0; h < vc.size(); h++){ //here we add the volt dependent current source
            int v1, v2, c1, c2;
            v1 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[0]);
            v2 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[1]);
            c1 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[2]);
            c2 = find_index(stand_posi_in_matrix, vc[h]->get_polarity()[3]);
            double trans = vc[h]->conductance(0).real();
            if(vc[h]->get_polarity()[0] != 0){
                if(vc[h]->get_polarity()[2] != 0){
                    mat(v1 + num_de, c1) -= trans;
                }
                if(vc[h]->get_polarity()[3] != 0){
                    mat(v1 +num_de, c2) += trans;
                }
            }
            if(vc[h]->get_polarity()[1] != 0){
                if(vc[h]->get_polarity()[2] != 0){
                    mat(v2 + num_de, c1) += trans;
                }
                if(vc[h]->get_polarity()[3] != 0){
                    mat(v2 + num_de, c2) -= trans;
                }
            }
        }
    }
    std::cout << "Mat" << std::endl;
    std::cout << mat << std::endl;
    Eigen::MatrixXd col_b (num, 1);
    for(int r = 0; r < num; r++){
        col_b (r, 0) = 0;
    }
    for(int r_n_volt = 0; r_n_volt < n_volt.size(); r_n_volt++){
        col_b (r_n_volt, 0) = v[n_volt[r_n_volt]]->get_volt().real();
    }
    for(int i = 0; i < a.size(); i++){ //to decide whether the current source is connected to the two nodes or one node connecting to the ground.
        if(a[i]->get_node()[0] == 0 && a[i]->dc_ac() == "dc"){
            if(a[i]->get_polarity()[0] != 0){
                col_b (find_index(stand_posi_in_matrix ,a[i]->get_polarity()[0]) + num_de, 0) = a[i]->get_current().real();
            }
            else{
                col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[1])+ num_de, 0) = - a[i]->get_current().real();
            }
        }
        if(a[i]->get_node()[0] != 0 && a[i]->dc_ac() == "dc"){
            col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[0])+ num_de, 0) = a[i]->get_current().real();
            col_b (find_index(stand_posi_in_matrix, a[i]->get_polarity()[1])+ num_de, 0) = - a[i]->get_current().real();
        }
    }
    for(int i = 0; i < ground_p.size(); i++){
        if(v[ground_p[i]]->get_node()[0] == 0){
            if(v[ground_p[i]]->get_polarity()[0] != 0){
                col_b (find_index(stand_posi_in_matrix, v[ground_p[i]]->get_node()[1]) + num_de, 0) = v[ground_p[i]]->get_volt().real();
            }
            if(v[ground_p[i]]->get_polarity()[1] != 0){
                col_b (find_index(stand_posi_in_matrix, v[ground_p[i]]->get_node()[1]) + num_de, 0) = - v[ground_p[i]]->get_volt().real();
            }
        }
    }
    std::cout << "c" << std::endl;
    std::cout << col_b << std::endl;
    Eigen::MatrixXd tmp_str = col_b;
    Eigen::MatrixXd guess_vot (num ,1);
    guess_vot = mat.colPivHouseholderQr().solve(col_b);
    Eigen::MatrixXd out (num - num_de, 1);
    for(int h = 0; h < out.rows(); h++){
        out(h, 0) = guess_vot(h,0);
    }
    return out;
}
bool compare_volt(Eigen::MatrixXd result, Eigen::MatrixXd test){
    bool decide = true;
    for(int i = 0; i < result.size(); i++){
        if(abs(result(i,0) - test(i,0)) > 0.0001){
           decide = false;
       }
    }
    return decide;
}
Eigen::MatrixXd get_standart_volt(std::vector<general*> input, int maxnode){
    Eigen::MatrixXd guess = build_guess_volt(input, maxnode);
    Eigen::MatrixXd iterate = build_iterate_matrix(input, maxnode, guess);
    Eigen::MatrixXd fvm = build_fvm_matrix(input, maxnode, guess);
    Eigen::MatrixXd result = guess;
    Eigen::MatrixXd testin (guess.rows(), 1);
    for(int g = 0 ; g< guess.rows(); g++){
        testin(g, 0) = 0;
    }
    int i = 0;
    std::vector<shortcircuit> s;
    std::vector<general*> test = short_circuit(input, s);
    while(!compare_volt(result , testin)){
        testin = result;
        iterate = build_iterate_matrix(input, maxnode, testin);
        fvm = build_fvm_matrix(input , maxnode, testin);
        result =testin - iterate.inverse() * fvm;   
    }
    result = recover_circuit(result,s );
    return result;
}

int main(){         
    std::vector <std::string> s;
    std::vector <general*> g;
    s = ReadInput("diode.txt");
    int maxn = 0;
    setting(s,g,maxn);
    
    Eigen::MatrixXd test;
    test = build_guess_volt(g, maxn);
    std::cout << test << std:: endl;
    Eigen::MatrixXd tmp;
    tmp = build_fvm_matrix(g, maxn, test);
    std::cout<< "Fvm" << std::endl;
    std::cout << tmp << std::endl;
    Eigen::MatrixXd final;
    final = get_standart_volt(g , maxn);
    std::cout << "final" << std::endl;
    std::cout << final << std::endl;
}
