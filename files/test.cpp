#include<iostream>
#include<complex>
#include<vector>
#include<cmath>
#include<string>
#include"part1.cpp"
#include<eigen3/Eigen/Dense>

void dc_volt (std::vector<voltsrc*> tmp, int &num, std::vector<int> &n, std::vector<int> &g){
    for (int i = 0; i < tmp.size(); i++){
            if(tmp[i]->get_node()[0] != 0 && tmp[i]->get_node()[1] != 0 && tmp[i]->dc_ac() == "dc"){
            num += 1;
            n.push_back(i);
            }
            if(tmp[i]->get_node()[0] == 0 && tmp[i]->dc_ac() == "dc"){
            g.push_back(tmp[i]->get_node()[1]);
            }
        }        
    }

void classify_comp(std::vector<general*> input, std::vector<Resistor*> &r, std::vector<bjt*> &b, std::vector<mosfet*> &m, std::vector<Diode*>&d, std::vector<voltsrc*> &v ){
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
    }
}

Eigen::MatrixXd build_guess_volt (std::vector<general*> input, int maxnode){
    std::vector<Resistor*> r;
    std::vector<bjt*> b;
    std::vector<mosfet*> m;
    std::vector<Diode*>d; 
    std::vector<voltsrc*>v;
    classify_comp(input, r, b, m, d, v);
     int num_de = 0;
    std::vector<int> g_volt;
// this vector is used to store the nodes of voltsrc connecting to ground.
    std::vector<int> n_volt;
// this vector contains the index for the voltsrc connecting to two nodes.
    dc_volt(v, num_de, n_volt, g_volt);
    int num = num_de + maxnode;
    Eigen::MatrixXd mat(num, num);
       for(int row = 0; row < num; row++){
        std::complex<double> tot_cond (0,0);
        std::vector<int > tmp2;
        tmp2.push_back(0);
        tmp2.push_back(row + 1);
        bool decide = false;       
        for(int c=0; c<num; c++){
            mat(row, c) = 0;
            }
        // construct model row full 0;
        for(int k = 0; k < g_volt.size(); k++){
        if(row == num_de + g_volt[k]-1){
            //connect to ground, put row
            // attention to the position of this.
            mat(g_volt[k]-1 + num_de, g_volt[k]-1) = 1;
            decide = true;
        }
        }
        if(row < num_de ){
            // connect to two nodes
            int posi = input[n_volt[row]]->get_polarity()[0];
            int nega = input[n_volt[row]]->get_polarity()[1];
            mat(row, posi-1) = 1;
            mat(row, nega-1) = -1;
            decide = true;                
            }
        std::complex<double> ground_cond (0,0);
        if(!decide){
        for(int index = 0; index < num; index++){
            std::vector<int> tmp;
            if(row - num_de < index){
                tmp.push_back((row-num_de+1));
                tmp.push_back((index+1));
            }
            else{
                tmp.push_back((index+1));
                tmp.push_back((row-num_de+1));
            }
            //construct compare vector;

            std::complex<double> cond (0,0);
            ground_cond = (0,0);
            for(int seq_comp = 0; seq_comp < input.size() ; seq_comp ++){
                if((input[seq_comp]-> nodes() == 2)&&(input[seq_comp]->get_node()== tmp)&&(input[seq_comp]->get_type() == "R")){
                    cond = cond - input[seq_comp]-> conductance(0);
                }
                //question : how to avoid this situation
                if((input[seq_comp]-> nodes() == 2) && (input[seq_comp]->get_node()==tmp2) && (input[seq_comp]->get_type() == "R")){
                    ground_cond = ground_cond + input[seq_comp]-> conductance(0);
                }
                //this function is used to find volt dependent current source.
                //translate the current output to the conductance at the two control node to the connect node.  
                if(input[seq_comp]->get_type() == "G") { 
                    if(input[seq_comp]->get_polarity()[0] == row){
                        if(input[seq_comp]->get_polarity()[2] == index){
                            cond = cond - input[seq_comp]-> conductance(0);
                        }
                        if(input[seq_comp]->get_polarity()[3] == index){
                            cond = cond + input[seq_comp]-> conductance(0);
                        }
                }
                    if(input[seq_comp]->get_polarity()[1] == row){
                        if(input[seq_comp]->get_polarity()[2] == index){
                            cond = cond + input[seq_comp]-> conductance(0);
                        }
                        if(input[seq_comp]->get_polarity()[3] == index){
                            cond = cond - input[seq_comp]-> conductance(0);
                        }
                    }
                    
                }
                

            }
            if(row - num_de != index){
                mat(row, index) = cond.real();
            }
            tot_cond = tot_cond - cond;
            }
            tot_cond += ground_cond;
            mat(row, row-num_de) = tot_cond.real();
            if(num_de > 0){
               for (int t = 0; t < num_de ; t++){
                   if(input[n_volt[t]]->get_polarity()[0] == (row - num_de + 1)){
                       mat(row, t + maxnode) = -1;
                   }
                   else if(input[n_volt[t]]->get_polarity()[1] == (row - num_de + 1)){
                       mat(row , t + maxnode) = 1;
                   }
               }
            }
        }
        }
    
    Eigen::MatrixXd col_b (num, 1);
    for(int r = 0; r < num; r++){
        col_b (r, 0) = 0;
    }
  

    for(int r_n_volt = 0; r_n_volt < n_volt.size(); r_n_volt++){
        col_b (r_n_volt, 0) = input[n_volt[r_n_volt]]->get_volt().real();
    }
    
    for(int i = 0; i < input.size(); i++){
        if(input[i]->get_type() == "I"){
            // to decide whether the current source is connected to the two nodes or one node connecting to the ground.
            if(input[i]->get_node()[0] == 0){
                if(input[i]->get_polarity()[0] != 0){
                    col_b (input[i]->get_polarity()[0] + num_de - 1, 0) = input[i]->get_current().real();
                }
                else{
                    col_b (input[i]->get_polarity()[1] + num_de - 1, 0) = - input[i]->get_current().real();
                }
            }
        else{
            col_b (input[i]->get_polarity()[0] + num_de -1, 0) = input[i]->get_current().real();
            col_b (input[i]->get_polarity()[1] + num_de -1, 0) = - input[i]->get_current().real();
        }
        }
        if((input[i]->get_type() == "V") && (input[i]->get_node()[0] == 0)){
            if(input[i]->get_polarity()[0] != 0){
                col_b (input[i]->get_node()[1] + num_de - 1, 0) = input[i]->get_volt().real();
            }
            if(input[i]->get_polarity()[1] != 0){
            col_b (input[i]->get_node()[1] + num_de - 1, 0) = - input[i]->get_volt().real();
            }
        }
    }
    Eigen::MatrixXd tmp_str = col_b;
            // for bjt    

    for(int i = 0; i < b.size(); i ++){
            // make assumption vbe = 0.7, vce =1;
            b[i]->initialize_bjt(0.7 , 0, 0.75);
            //for node C
            col_b(b[i]->get_polarity()[0]+ num_de -1,0) -= b[i]->ic();
            // for node B
            col_b(b[i]->get_polarity()[1] + num_de -1,0) -= b[i]->ib();
            // for node E
            col_b(b[i]->get_polarity()[2] + num_de -1,0) += b[i]->ie();
    }
            // for mosfet
        for(int i = 0; i < m.size(); i++){
            //make assumption vgs = vt, vds = 0;
            m[i]->initialize_mos(m[i]->vt(), 0, 0);
            //for node D 
            col_b(m[i]->get_polarity()[2] + num_de -1, 0) -= m[i]->id();
            //for node S
            col_b(m[i]->get_polarity()[1] + num_de -1, 0) += m[i]->id();
        }

        for(int n = 0; n < d.size(); n++){
            col_b(d[n]->get_polarity()[0] + num_de -1, 0) -= d[n]->id(0.7, 0);
            col_b(d[n]->get_polarity()[1] + num_de -1, 0) += d[n]->id(0.7, 0);
        }

    Eigen::MatrixXd guess_vot (num ,1);
    guess_vot = mat.colPivHouseholderQr().solve(col_b);
    bool check_precise = true;
    for (int h = 0; h < b.size(); h++){
        if((guess_vot(b[h]->get_polarity()[1] + num_de -1, 0) - guess_vot(b[h]->get_polarity()[2] + num_de -1, 0) > 1 )|| (guess_vot(b[h]->get_polarity()[1] + num_de -1 , 0) - guess_vot(b[h]->get_polarity()[2] + num_de -1, 0) < 0)){
            check_precise = false;
        }
    }
    if(!check_precise){
        guess_vot = mat.colPivHouseholderQr().solve(tmp_str);
    }
    
    return guess_vot;
}

Eigen::MatrixXd build_iterate_matrix (std::vector<general*> input, int maxnode, std::vector<double> volt_vector){

std::vector<Resistor*> r;
std::vector<bjt*> b;
std::vector<mosfet*> m;
std::vector<Diode*>d; 
std::vector<voltsrc*>v;
classify_comp(input, r, b, m, d, v);
int num_de = 0;
    std::vector<int> g_volt;
// this vector is used to store the nodes of voltsrc connecting to ground.
    std::vector<int> n_volt;
// this vector contains the index for the voltsrc connecting to two nodes.
dc_volt(v, num_de, n_volt, g_volt);
    int num = num_de + maxnode;
    Eigen::MatrixXd mat(num, num);


    for(int row = 0; row < num; row++){
       double tot_cond= 0;
        std::vector<int > tmp2;
        tmp2.push_back(0);
        tmp2.push_back(row + 1);
        bool decide = false;
       for(int c=0; c<num; c++){
            mat(row, c) = 0;
            }
        // construct model row full 0;
        for(int k = 0; k < g_volt.size(); k++){
        if(row == num_de + g_volt[k]-1){
            //connect to ground, put row
            // attention to the position of this.
            mat(g_volt[k]-1 + num_de, g_volt[k]-1) = 1;
            decide = true;
        }
        }
        if(row < num_de ){
            // connect to two nodes
            int posi = input[n_volt[row]]->get_polarity()[0];
            int nega = input[n_volt[row]]->get_polarity()[1];
            mat(row, posi-1) = 1;
            mat(row, nega-1) = -1;
            decide = true;                
            }
        double ground_cond= 0;
        if(!decide){
        for(int index = 0; index < num; index++){
            std::vector<int> tmp;
            if(row - num_de < index){
                tmp.push_back((row-num_de+1));
                tmp.push_back((index+1));
            }
            else{
                tmp.push_back((index+1));
                tmp.push_back((row-num_de+1));
            }
            //construct compare vector;

            double cond = 0;
            ground_cond =0;
            for(int seq_comp = 0; seq_comp < input.size() ; seq_comp ++){
                if((input[seq_comp]->get_node()== tmp)&&(input[seq_comp]->get_type() == "R")){
                    cond = cond - input[seq_comp]-> conductance(0).real();
                }
                //question : how to avoid this situation
                if((input[seq_comp]->get_node()==tmp2) && (input[seq_comp]->get_type() == "R")){
                    ground_cond = ground_cond + input[seq_comp]-> conductance(0).real();
                }
                //this function is used to find volt dependent current source.
                //translate the current output to the conductance at the two control node to the connect node.
                if(input[seq_comp]->get_type() == "L" && input[seq_comp]->get_node() == tmp){
                    cond = cond - pow(10, 13);
                }  
                if(input[seq_comp]->get_type() == "G") { 
                    if(input[seq_comp]->get_polarity()[0] == row){
                        if(input[seq_comp]->get_polarity()[2] == index){
                            cond = cond - input[seq_comp]-> conductance(0).real();
                        }
                        if(input[seq_comp]->get_polarity()[3] == index){
                            cond = cond + input[seq_comp]-> conductance(0).real();
                        }
                }
                    if(input[seq_comp]->get_polarity()[1] == row){
                        if(input[seq_comp]->get_polarity()[2] == index){
                            cond = cond + input[seq_comp]-> conductance(0).real();
                        }
                        if(input[seq_comp]->get_polarity()[3] == index){
                            cond = cond - input[seq_comp]-> conductance(0).real();
                        }
                    }
                    
                }
                

            }
            if(row - num_de != index){
                mat(row, index) = cond;
            }
            tot_cond = tot_cond - cond;
            }
            tot_cond += ground_cond;
            mat(row, row-num_de) = tot_cond;
            if(num_de > 0){
               for (int t = 0; t < num_de ; t++){
                   if(input[n_volt[t]]->get_polarity()[0] == (row - num_de + 1)){
                       mat(row, t + maxnode) = -1;
                   }
                   else if(input[n_volt[t]]->get_polarity()[1] == (row - num_de + 1)){
                       mat(row , t + maxnode) = 1;
                   }
               }
            }
        }
        }

        // for diode
        for (int i = 0; i < d.size(); i++)
        {   //anode row
            mat(d[i]->get_polarity()[0] -1 + num_de, d[i]->get_polarity()[1] -1) += - d[i]->g(volt_vector[d[i]->get_polarity()[0]], volt_vector[d[i]->get_polarity()[1]]);
            mat(d[i]->get_polarity()[0] -1 + num_de, d[i]->get_polarity()[0] -1) += d[i]->g(volt_vector[d[i]->get_polarity()[0]], volt_vector[d[i]->get_polarity()[1]]);
            //cathode row
            mat(d[i]->get_polarity()[1] -1 + num_de, d[i]->get_polarity()[1] -1) += d[i]->g(volt_vector[d[i]->get_polarity()[0]], volt_vector[d[i]->get_polarity()[1]]);
            mat(d[i]->get_polarity()[1] -1 + num_de, d[i]->get_polarity()[0] -1) += -d[i]->g(volt_vector[d[i]->get_polarity()[0]], volt_vector[d[i]->get_polarity()[1]]);
        }

        // for bjt

        for(int m = 0; m < b.size(); m++){
            b[m]->initialize_bjt(volt_vector[b[m]->get_polarity()[0]], volt_vector[b[m]->get_polarity()[1]], volt_vector[b[m]->get_polarity()[2]]);

            // for bjt with b on that row.
            mat(b[m]->get_polarity()[1] -1 + num_de,b[m]->get_polarity()[1] - 1) += b[m]->b_b(); 
            mat(b[m]->get_polarity()[1] -1 + num_de, b[m]->get_polarity()[2] - 1) += b[m]->b_e();
            mat(b[m]->get_polarity()[1] -1 + num_de, b[m]->get_polarity()[0] -1) += b[m]->b_c();
            
            // for bjt with e on that row.
            mat(b[m]->get_polarity()[2] -1, b[m]->get_polarity()[2]-1) += b[m]->e_e(); 
            mat(b[m]->get_polarity()[2] -1, b[m]->get_polarity()[1]-1) += b[m]->e_b();
            mat(b[m]->get_polarity()[2] -1, b[m]->get_polarity()[0]-1) += b[m]->e_c();
            
            // for bjt with c on that row
            mat(b[m]->get_polarity()[0] -1, b[m]->get_polarity()[2]-1) -= b[m]->c_e(); 
            mat(b[m]->get_polarity()[0] -1, b[m]->get_polarity()[1]-1) += b[m]->c_b();
            mat(b[m]->get_polarity()[0] -1, b[m]->get_polarity()[0]-1) += b[m]->c_c();
        }
        // for mosfet

        for(int g = 0; g < m.size(); g++){
            m[g]->initialize_mos(volt_vector[m[g]->get_polarity()[0]], volt_vector[m[g]->get_polarity()[1]], volt_vector[m[g]->get_polarity()[2]]);
            
            // for mos with s on that row
            
            mat(m[g]->get_polarity()[1] -1, m[g]->get_polarity()[0]-1) += m[g]->Id_g();
            mat(m[g]->get_polarity()[1] -1, m[g]->get_polarity()[1]-1) += m[g]->Id_s();
            mat(m[g]->get_polarity()[1] -1, m[g]->get_polarity()[2]-1) += m[g]->Id_d();

            // for mos with d on that row
            mat(m[g]->get_polarity()[2]-1, m[g]->get_polarity()[0] -1) -= m[g]->Id_g();
            mat(m[g]->get_polarity()[2]-1, m[g]->get_polarity()[1] -1) -= m[g]->Id_g();
            mat(m[g]->get_polarity()[2]-1, m[g]->get_polarity()[2] -1) -= m[g]->Id_g();
        }

        //

        return mat;

}


/*

int main(){         
    vector <string> s;
    vector <general*> g;
    s = ReadInput("input3.txt");
    int maxn = 0;
    setting(s,g,maxn);
    
    Eigen::MatrixXd guess = build_guess_volt(g, maxn);
    std::vector<double> tmp;
    for(int i = 0; i < guess.rows(); i++){
        tmp.push_back(guess(i, 0));
    }

    Eigen::MatrixXd test = build_iterate_matrix(g, maxn, tmp);
    std::cout << test << std::endl;
    std::cout << guess << std::endl;
}
*/