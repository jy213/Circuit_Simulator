#include<iostream>
#include "part1.cpp"
#include<complex>
#include <string>
#include <eigen3/Eigen/Dense>

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

void output(vector<double> fre, vector<double> sol){
    ofstream outFile;
    outFile.open("output.txt");
    for(int i=0;i<fre.size();i++){
      outFile << fre[i] << "," << sol[i] << endl;  
    }
    outFile.close();
}

void num_of_devolt (std::vector<general*> input, int &num, std::vector<int> &n, std::vector<int> &g){
    for (int i = 0; i < input.size(); i++){
        if((input[i]->get_type() == "V") && (input[i]->get_node()[0] != 0 && input[i]->get_node()[1] != 0)){
            num += 1;
            n.push_back(i);
        }
        if((input[i]->get_type() == "V") && (input[i]->get_node()[0] == 0)){
            g.push_back(input[i]->get_node()[1]);
        }
    }
}
// To build conductance Matrix
Eigen::MatrixXcd build_cond_matrix(std::vector<general*> input, int maxnode, double frequency){
    int num_de = 0;
    std::vector<int> g_volt;
// this vector is used to store the nodes of voltsrc connecting to ground.
    std::vector<int> n_volt;
// this vector contains the index for the voltsrc connecting to two nodes.
 num_of_devolt(input, num_de, n_volt, g_volt);
    int num = num_de + maxnode;
    Eigen::MatrixXcd mat(num, num);
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
                if((input[seq_comp]-> nodes() == 2)&&(input[seq_comp]->get_node()== tmp)){
                    cond = cond - input[seq_comp]-> conductance(frequency * 2 * PI);
                }
                //question : how to avoid this situation
                if((input[seq_comp]-> nodes() == 2) && (input[seq_comp]->get_node()==tmp2)){
                    ground_cond = ground_cond + input[seq_comp]-> conductance(frequency * 2 * PI);
                }
                //this function is used to find volt dependent current source.
                //translate the current output to the conductance at the two control node to the connect node.  
                if(input[seq_comp]->get_type() == "G") { 
                    if(input[seq_comp]->get_polarity()[0] == row){
                        if(input[seq_comp]->get_polarity()[2] == index){
                            cond = cond - input[seq_comp]-> conductance(frequency * 2 * PI);
                        }
                        if(input[seq_comp]->get_polarity()[3] == index){
                            cond = cond + input[seq_comp]-> conductance(frequency * 2 * PI);
                        }
                }
                    if(input[seq_comp]->get_polarity()[1] == row){
                        if(input[seq_comp]->get_polarity()[2] == index){
                            cond = cond + input[seq_comp]-> conductance(frequency * 2 * PI);
                        }
                        if(input[seq_comp]->get_polarity()[3] == index){
                            cond = cond - input[seq_comp]-> conductance(frequency * 2 * PI);
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
    // the end of the row loop;
    // to build the right col with 1 , -1.
        }
    // the end of the matrix loop;
    return mat;
}
Eigen::MatrixXcd build_b (std::vector<general*> input, int maxnode){
        int num_de = 0;
    std::vector<int> g_volt;
// this vector is used to store the nodes of voltsrc connecting to ground.
    std::vector<int> n_volt;
// this vector contains the index for the voltsrc connecting to two nodes.
 num_of_devolt(input, num_de, n_volt, g_volt);
    int num = num_de + maxnode;
    Eigen::MatrixXcd col_b (num, 1);
    for(int r = 0; r < num; r++){
        col_b (r, 0) = 0;
    }

    for(int r_n_volt = 0; r_n_volt < n_volt.size(); r_n_volt++){
        col_b (r_n_volt, 0) = input[n_volt[r_n_volt]]->get_volt();
    }
    
    for(int i = 0; i < input.size(); i++){
        if(input[i]->get_type() == "I"){
            // to decide whether the current source is connected to the two nodes or one node connecting to the ground.
            if(input[i]->get_node()[0] == 0){
                if(input[i]->get_polarity()[0] != 0){
                    col_b (input[i]->get_polarity()[0] + num_de - 1, 0) = input[i]->get_current();
                }
                else{
                    col_b (input[i]->get_polarity()[1] + num_de - 1, 0) = - input[i]->get_current();
                }
        }
        else{
            col_b (input[i]->get_polarity()[0] + num_de -1, 0) = input[i]->get_current();
            col_b (input[i]->get_polarity()[1] + num_de -1, 0) = - input[i]->get_current();
        }
        }
        if((input[i]->get_type() == "V") && (input[i]->get_node()[0] == 0)){
            if(input[i]->get_polarity()[0] != 0){
                col_b (input[i]->get_node()[1] + num_de - 1, 0) = input[i]->get_volt();
            }
            if(input[i]->get_polarity()[1] != 0){
            col_b (input[i]->get_node()[1] + num_de - 1, 0) = - input[i]->get_volt();
            }
        }
    }

    return col_b;
}

Eigen::VectorXcd find_solution(Eigen::MatrixXcd gmatrix, Eigen::VectorXcd vicolumn) {
	Eigen::VectorXcd Vx = gmatrix.colPivHouseholderQr().solve(vicolumn);
	return Vx;
}



int main(){
    vector <string> s;
    vector <general*> g;
    s = ReadInput("input2.txt");
    int maxn = 0;
    setting(s,g,maxn);
    Eigen::MatrixXcd test;
    test = build_b(g, maxn);
    std::cout << test << std::endl;
    Eigen::MatrixXcd test2;
    test2 = build_cond_matrix(g, maxn, 100);
    std::cout << test2 << std::endl;
    Eigen::MatrixXcd sol;
    sol = find_solution(test2, test);
    std::cout << sol<< std::endl;
}


        