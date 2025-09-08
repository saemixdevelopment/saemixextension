#include <Rcpp.h>
#include <unordered_set>
#include <vector>
#include <string> 
#include <cmath>
using namespace Rcpp; 

std::vector<float> import_rgrid_b(const NumericVector& mat , IntegerVector& dims , std::vector<std::pair<int,int>>& vec_pad , std::vector<int>& padded_strides,
    int d, std::vector<int>& padded_dims  , float fill = 0.0){

    std::vector<float> result ; 
    int total_size = 1 ; // permettra de calculer la taille du vecteur 
    std::vector<int> strides(d) ; 
    strides[0] = 1 ; 
    padded_strides[0] = 1 ; 

    for (int i = 0 ; i < d ; ++i){
        padded_dims[i] = dims[i] + vec_pad[i].first + vec_pad[i].second ; // Pour chaque sous vecteur , on met sa dimension original + l'allocation des deux côtés
        total_size*= padded_dims[i] ;  
        if (i > 0){
            strides[i] = strides[i-1] * dims[i - 1] ; 
            padded_strides[i] = padded_strides[i-1] * padded_dims[i-1] ; 
        }
    }

    if (total_size == mat.size()){
        result.assign(mat.begin() , mat.end());
    } else {
        result.resize(total_size,fill) ; // on fait notre nouveau vecteur avec la bonne allocation. 
        std::vector<int> index(d) ; // Base qui nous serviras à l'avenir pour marquer les index 

        for(int lin_index = 0 ; lin_index < mat.size() ; ++lin_index){ // pour chaque index de notre matrice original 
            int remaining = lin_index ; // reste des opérations pour retrouver à quelle position nous sommes
            for(int i = d - 1 ; i >=0 ; --i){ // pour chaque dimension 
                index[i] = remaining / strides[i] ; // on récupère le reste car index est un int , numéro de la dimension où l'on est 
                remaining %= strides[i] ; // dans quelle position dans la dimension d - 1  , nous sommes 
            }
            int padded_index = 0 ; 
            for (int i = 0 ; i < d ; ++i){
                int shifted = index[i] + vec_pad[i].first ; // Permet de shift en prenant en compte le pad de la dimension actuel 
                padded_index += shifted * padded_strides[i] ; // fait le cumule des index_pad pour retrouver sa position dans le vecteur liénaire. 
            }
            result[padded_index] = mat[lin_index] ; // On place aux bonnes positions les index que l'on avait entre la matrice pad et non pad. 
        }   
    }
    return result ;  
}


std::vector<std::vector<std::vector<float>>> import_rgrid_x(List& input_list , std::vector<std::vector<std::pair<int,int>>>& vec_pad , float fill = 0.0){
    int global_size = input_list.size() ; 
    std::vector<std::vector<std::vector<float>>> output(global_size) ; 
   
    for (int i = 0 ;  i < global_size ; i++){
        const List& sub_list = input_list[i] ; 
        std::vector<std::vector<float>> padded_inner(sub_list.size()) ; 

        for (int j = 0 ; j < sub_list.size() ; j++ ){
            const NumericVector& vec = sub_list[j] ; 
            int original_vec_size = vec.size() ; 
            int total_size = original_vec_size + vec_pad[i][j].first + vec_pad[i][j].second ; 
            std::vector<float> padded_vec(total_size,fill) ; 
            for (int k = 0 ; k < original_vec_size ; ++k){
                padded_vec[vec_pad[i][j].first + k] = vec[k] ; 
            }
            padded_inner[j] = padded_vec ; 
        }
     output[i] = padded_inner ; 
    }
    return output ; 
}

std::vector<int> get_unique(std::vector<int>& vec){
    std::unordered_set<int> set ; 
    std::vector<int> vec_output ; 
    for (int i = 0 ; i < vec.size() ; i++ ){
        if (set.insert(vec[i]).second){ 
            vec_output.push_back(vec[i]) ; 
        }
    }
    return (vec_output) ; 
}

int index_detect(std::vector<float>& vec,float& limit){
    for (int i = 0 ; i < vec.size() ; i++){
        if (vec[i] > limit){
            return (i) ; 
        }
    }
    return -1  ; 
}

void expand_vector(std::vector<std::pair<int,int>>& vec_combi , std::vector<int>& combi
     , std::vector<std::vector<int>>& vec_output , int depth = 0){

        if (depth == vec_combi.size()){
            vec_output.push_back(combi) ; 
            return ; 
        }

        combi[depth] = vec_combi[depth].first ; 
        expand_vector(vec_combi , combi , vec_output , depth + 1) ; 

        combi[depth] = vec_combi[depth].second ; 
        expand_vector(vec_combi , combi , vec_output , depth+1 ) ; 

}

void expand_vector_float(std::vector<std::pair<float,float>>& vec_combi , std::vector<float>& combi
    , std::vector<std::vector<float>>& vec_output , int depth = 0){

       if (depth == (vec_combi.size())){
           vec_output.push_back(combi) ; 
           return ; 
       }

       combi[depth] = vec_combi[depth].first ; 
       expand_vector_float(vec_combi , combi , vec_output , depth + 1) ; 

       combi[depth] = vec_combi[depth].second ; 
       expand_vector_float(vec_combi , combi , vec_output , depth +1 ) ; 

}

int get_index(std::vector<int>& dims , std::vector<int>& index , std::vector<std::vector<std::pair<int,int>>>& vec_pad ,
    std::vector<int>& j_padded_strides , std::vector<std::vector<std::pair<int,int>>> vec_pad_left_right , int& j ){
    int d = dims.size() ; 
    int lin_index = 0 ; 
    for (int i = 0 ; i < d ; ++i){
        lin_index+= (index[i] + vec_pad[j][i].first - vec_pad_left_right[j][i].first) * j_padded_strides[i] ; 
        // + 1 ici car les index sur C++ commencent par zéros , or , comme on multiplie , on veut la somme des dimensions
    }
    return lin_index ; 
}

void get_borne_min_max(std::vector<std::vector<std::pair<int,int>>>& vec_pad, NumericMatrix& psi , 
    std::vector<int>& unique_id , int& len_rgrid_x , List& mat_rgrid_x , int& d , int& N , double& alg , 
    std::vector<int>& transfor_parj){

    int count = 1 ; 
    float value ; 
    float U0 ; 
    float U1 ; 
    int j ; 

    for (int k = 0 ; k < psi.nrow() ; k++){
        j = ((unique_id[k]-1) % len_rgrid_x) ;

        const List& sublist_j = mat_rgrid_x[j] ;  
        const NumericVector& psi_k = psi.row(k) ; 
        for (int l = 0 ; l < d ; l++){
            const NumericVector& vec_l = sublist_j[l];
            if (psi_k[l] < vec_l[0]){ 
                U0 = vec_l[0] ; 
                U1 = vec_l[1]  ;
                if (transfor_parj[l] == 0){
                    value = U0 - count * (U1 - U0) ; 
                    while (psi_k[l] < value){
                        count += 1 ;
                        value = U0 - count * (U1 - U0) ;       
                    }
                } else {
                    value = std::pow(U0,count+1)/std::pow(U1,count) ; 
                    while (psi_k[l] < value){
                        count += 1 ; 
                        value = std::pow(U0,count+1)/std::pow(U1,count) ; 
                    }
                }
                if (vec_pad[j][l].first < count){
                    vec_pad[j][l].first = count ; 
                }
                count = 1 ; 
                continue ; 
            }
            if (psi_k[l] > vec_l[vec_l.size() - 1]){ 
                U0 = vec_l[vec_l.size() - 1] ; 
                U1 = vec_l[vec_l.size() - 2]  ;
                if (transfor_parj[l] == 0){
                    value = U0 + count * (U0 - U1) ;  
                    while (psi_k[l] > value){
                        count += 1 ; 
                        value = U0 + count * (U0 - U1) ; 
                    }
                } else {
                    value =  std::pow(std::pow(U0,alg) + count * (std::pow(U0,alg) - std::pow(U1,alg)),1/alg) ; 
                    while (psi_k[l] > value){ 
                        count += 1 ; 
                        value = std::pow(std::pow(U0,alg) + count * (std::pow(U0,alg) - std::pow(U1,alg)),1/alg) ;     
                    }
                }
                if (vec_pad[j][l].second < count){
                    vec_pad[j][l].second = count ; 
                }
                count = 1 ; 
                continue ; 
            }
        }
    }
}

void get_size_each_id(std::vector<int>& vec_size_each_id ,std::vector<int>& id , std::vector<int>&strides_id){
    int count = 1 ; 
    int j = 0 ; 
    int strides = 0  ; 
    strides_id[0] = 0 ; 
    for (int i = 0 ; i < (id.size() - 1) ; i++){
        if (id[i] == id[i+1] && i < ((id.size() - 2))){
            count+= 1 ; 
        } else {
            if (i == ((id.size() - 2))){
                count+=1 ; 
            }
            vec_size_each_id[j] = count ; 
            strides += count ; 
            strides_id[j+1] = strides ; 
            count = 1 ; 
            j += 1 ; 
        }
    }
}

std::vector<float> produit_matriciel_poids_sim(std::vector<std::vector<float>>& A ,std::vector<float>& B ){
    // La matrice B est une matrice qui contient forcemment n rows et 1 colonne 
    int nrow_a = A.size() ; 
    int  ncol_a = A[0].size() ; 
    
    int nrow_b = B.size() ; 
    int ncol_b = 1 ; 

    int size_new_matrix = ncol_a * ncol_b ; 
    std::vector<float> output_matrix(size_new_matrix) ; 

    float valeur = 0 ; 
    for (int row_mat = 0 ; row_mat < size_new_matrix ; row_mat++){
        for (int row_a = 0 ; row_a < nrow_a ; row_a++){
            valeur += A[row_a][row_mat] * B[row_a] ;   
        }
        output_matrix[row_mat] = valeur ; 
        valeur = 0;  
    }
    return output_matrix ; 
}

// [[Rcpp::export]]
List main_function(List mat_rgrid_b , List mat_rgrid_x , List mat_rgrid_f , NumericMatrix psi, 
std::vector<int> id , Function compute_simulation , std::vector<int> transfor_parj){

  
    int len_rgrid_x = mat_rgrid_x.size() ;  
    int N = psi.nrow() ; 
    int d = psi.ncol() ;
    std::vector<std::vector<std::pair<int,int>>> vec_pad_left_right(len_rgrid_x , 
        std::vector<std::pair<int,int>>(d,std::make_pair(0,0))) ; 
    std::vector<int> unique_id = get_unique(id) ; 
    std::vector<std::vector<std::pair<int,int>>> vec_pad(len_rgrid_x , 
        std::vector<std::pair<int,int>>(d,std::make_pair(0,0))) ; 
    double alg = 0.5 ;  
  
    get_borne_min_max(vec_pad,psi,unique_id,len_rgrid_x,mat_rgrid_x,d,N,alg,transfor_parj) ; 
     
    std::vector<std::vector<float>> pad_mat_rgrid_b ; 
    std::vector<std::vector<int>> list_padded_strides ; 
    std::vector<std::vector<int>> list_padded_dims ; 
    std::vector<int> strides_id(unique_id.size()+1) ; 
    std::vector<int> vec_size_each_id(unique_id.size()) ; 
    get_size_each_id(vec_size_each_id,id,strides_id) ; 

    for (int i = 0 ; i < mat_rgrid_b.size() ; i++){
        const NumericVector& mat = mat_rgrid_b[i] ; 
        IntegerVector dims = mat.attr("dim") ; 
        std::vector<int> padded_strides(dims.size()) ; 
        int d = dims.size() ; // on récupère le nombre de dimension qu'on a dans notre jeu de donnée
        std::vector<int> padded_dims(d) ;  // on fait un vecteur contenant 5 sous Vecteur
        std::vector<float> padded_mat = import_rgrid_b(mat,dims,vec_pad[i],padded_strides,d,padded_dims) ; 
        pad_mat_rgrid_b.push_back(padded_mat) ; 
        list_padded_strides.push_back(padded_strides) ; 
        list_padded_dims.push_back(padded_dims) ;         
    }

  
    std::vector<std::vector<std::vector<float>>> pad_rgrid_x = import_rgrid_x(mat_rgrid_x,vec_pad) ; 
    std::vector<std::vector<float>> vec_psi_a_visiter(len_rgrid_x) ;
    std::vector<std::vector<float>> vec_poids_id(N) ; 
    std::vector<std::vector<int>> vec_index_a_visiter(N) ; 
    bool new_value = false ; 

    for (int k  = 0 ; k < N ; k++){
        int j = ((unique_id[k]-1) % len_rgrid_x) ; 

        std::vector<float> psi_k(psi.row(k).begin() , psi.row(k).end()) ; 
        std::vector<std::vector<float>>& j_rgrid_x = pad_rgrid_x[j] ; 
        std::vector<float>& j_rgrid_b = pad_mat_rgrid_b[j] ; 
        std::vector<int>& j_padded_strides = list_padded_strides[j] ; 
        std::vector<int>& j_dims = list_padded_dims[j] ; 

        std::vector<std::pair<int,int>> ii(d) ; 
        std::vector<std::pair<float,float>> pii(d) ; 

        float pil ; 

        for (int l = 0 ; l < d ; l++){
            std::vector<float>& jl_rgrid_x = j_rgrid_x[l] ; 
            int& left = vec_pad_left_right[j][l].first ;
            int& right = vec_pad_left_right[j][l].second ; 
            int i2l = index_detect(jl_rgrid_x , psi_k[l]) ; 
            if (i2l > -1){
                i2l = i2l - vec_pad[j][l].first + left  ; 
            }
            if (i2l == -1){
                int ni = jl_rgrid_x.size() - 1 ; 
                int n_max  =  ni - vec_pad[j][l].second + right  ;
                while(psi_k[l] > jl_rgrid_x[n_max]){ 
                    if (transfor_parj[l] == 0){
                        jl_rgrid_x[n_max+1] = jl_rgrid_x[n_max] + (jl_rgrid_x[n_max] - jl_rgrid_x[n_max-1]) ; 
                    } else {
                        jl_rgrid_x[n_max+1] = std::pow(2*(std::pow(jl_rgrid_x[n_max],alg)) - (std::pow(jl_rgrid_x[n_max-1],alg)),1/alg) ; 
                    }
                    right+=1 ; 
                    n_max+=1 ; 
                }
                i2l = ni - (vec_pad[j][l].first - left) - (vec_pad[j][l].second - right)  ; 
            }
            if (i2l == 0){
                while (psi_k[l] < jl_rgrid_x[vec_pad[j][l].first-left]){
                    if (transfor_parj[l] == 0){
                        jl_rgrid_x[vec_pad[j][l].first-1-left] = jl_rgrid_x[vec_pad[j][l].first-left] + (jl_rgrid_x[vec_pad[j][l].first-left] - jl_rgrid_x[vec_pad[j][l].first-left+1]) ; 
                    } else {
                        jl_rgrid_x[vec_pad[j][l].first-1-left] = std::pow(jl_rgrid_x[vec_pad[j][l].first-left],2)/jl_rgrid_x[vec_pad[j][l].first-left+1] ; 
                    }
                    left+=1 ; 
                }
                i2l = 1 ; 
            }
            ii[l] = std::pair<int,int>(i2l-1,i2l) ; 
            pil = (jl_rgrid_x[i2l+vec_pad[j][l].first-left] - psi_k[l]) / (jl_rgrid_x[i2l+vec_pad[j][l].first-left] - jl_rgrid_x[i2l-1+vec_pad[j][l].first-left]) ; 
            pii[l] = std::pair<float,float>(pil,1-pil) ; 
        }
        std::vector<std::vector<int>> ig ; 
        std::vector<int> current(d) ; 
        expand_vector(ii  , current , ig) ; 
        std::vector<std::vector<float>> expand_whi ; 
        std::vector<float> current_float(d) ; 
        expand_vector_float(pii , current_float, expand_whi) ; 
        float prod ; 
        std::vector<float> whi ;
        // Permet simplement de multiplier tous les élements de la même ligne d'une matrice 
        for (int i = 0 ; i < expand_whi.size() ; i++){
            prod = 1 ; 
            for (int l = 0 ; l < expand_whi[i].size() ; l++ ){
                prod = prod * expand_whi[i][l] ; 
            }
            whi.push_back(prod) ; 
        }

        vec_poids_id[k] = whi ; 

        for (int l = 0 ; l < ig.size() ; l++){
            std::vector<int> igl = ig[l] ; 
            int index = get_index(j_dims,igl,vec_pad,j_padded_strides,vec_pad_left_right,j) ;  
            for (int p = 0 ; p < d ; p++){
                igl[p] = igl[p] + vec_pad[j][p].first -  vec_pad_left_right[j][p].first ; 
            } 
            int value = j_rgrid_b[index] ; 
            if (value == 0){ 
                std::vector<float> psi_m ; 
                psi_m.push_back(j_rgrid_x[0][igl[0]]) ; 
                if (d >= 2) {
                    for (int m = 1 ; m < d ; m++){ 
                        psi_m.push_back(j_rgrid_x[m][igl[m]]) ; 
                    }
                }
                vec_psi_a_visiter[j].insert(vec_psi_a_visiter[j].end() , psi_m.begin() , psi_m.end());
                vec_index_a_visiter[k].push_back(index) ; 
                // permettra de switch sur pksim et de faire les simulations
                new_value = true ; 
                j_rgrid_b[index] = -1 ; 
            } else {
                vec_index_a_visiter[k].push_back(index)  ; 
            }
        }
    }
    
   
    List input ; 

    if (new_value){
         input =  compute_simulation(wrap(vec_psi_a_visiter),id) ;  // permet d'intéragir avec R et récupérer les résultats des simulations à partir des psi. 
    }
    // deuxième partie avec récupération des résultats. 

    int nb_chains = N / len_rgrid_x  ; 
    int k ; 
    int total_strides = strides_id[unique_id.size()] ;  // pas besoin de faire - 1 car on a rajouté un zéro à la position 0
    NumericMatrix* j_rgrid_f = nullptr ; 
    NumericMatrix* pre_alloc_grid = nullptr;
    std::vector<float> Sim_prediction(total_strides) ;

    int ncol_f ; 
    int nrow_f ; 

    for (int j = 0 ; j < len_rgrid_x ; j++){
        int pre_alloc_for_j = 0 ; 

        if (!Rf_isNull(mat_rgrid_f[j]) && Rf_isMatrix(mat_rgrid_f[j])){
            j_rgrid_f = new NumericMatrix(mat_rgrid_f[j]) ; 
            IntegerVector j_dim = (*j_rgrid_f).attr("dim") ; 
            ncol_f =  j_dim[1] + 1 ; 
        } else {
            ncol_f = 1 ;  // je mets à 1 , pour que ça soit comme en R 
        }
        nrow_f = vec_size_each_id[j] ; 
        if (new_value && !vec_psi_a_visiter[j].empty()){
            pre_alloc_for_j = vec_psi_a_visiter[j].size() / d ; 
            pre_alloc_grid = new NumericMatrix(nrow_f,ncol_f+pre_alloc_for_j-1) ; 
    
            for (int row = 0 ; row < nrow_f ; row++){
                for (int col = 0 ; col < (ncol_f-1) ; col++){
                    (*pre_alloc_grid)(row,col) = (*j_rgrid_f)(row,col) ; 
                }
            }
            delete j_rgrid_f ; 
            j_rgrid_f = pre_alloc_grid ;  
        }
        std::vector<float>& j_rgrid_b = pad_mat_rgrid_b[j] ; 
        int iterator_j = 0 ; 

        for (int chains = 0 ; chains < nb_chains ; chains++){ 
            k = j  + len_rgrid_x * chains ;    
            int iterator_index = 0 ; 
            std::vector<float>& whi_k = vec_poids_id[k] ;  
            std::vector<std::vector<float>> sim_Results_k(vec_index_a_visiter[k].size()) ; 

            for (const int& index : vec_index_a_visiter[k]){
                if (j_rgrid_b[index] > 0 ){
                        NumericVector col_vec = (*j_rgrid_f)(_,(j_rgrid_b[index]-1)) ; 
                        sim_Results_k[iterator_index] = as<std::vector<float>>(col_vec); 
                } else {
                    NumericVector vec_result = as<NumericVector>(input[j])[Range(0 + iterator_j * vec_size_each_id[k],vec_size_each_id[k] + vec_size_each_id[k] * iterator_j - 1)] ;
                    iterator_j+= 1 ;   
                    sim_Results_k[iterator_index] = as<std::vector<float>>(vec_result) ; 
                    for (int i = 0; i < (*j_rgrid_f).nrow(); ++i) {
                        (*j_rgrid_f)(i, (ncol_f-1)) = vec_result[i];
                    }
                    j_rgrid_b[index] = static_cast<float>(ncol_f) ; 
                    ncol_f += 1 ; 
                   // il faudra à l'avenir gérer quand la valeur du résultat a eu une erreur.
                }
                iterator_index+=1 ; 
            }
            std::vector<float> result = produit_matriciel_poids_sim(sim_Results_k,whi_k) ; 
            int i = 0 ; 
            for (int m =  strides_id[k] ; m < strides_id[k+1] ; m++){
                Sim_prediction[m] = result[i] ; 
                i+= 1  ; 
            }
        }
        if(pre_alloc_grid != nullptr){
            mat_rgrid_f[j] = (*j_rgrid_f) ; 
            delete pre_alloc_grid ; 
            pre_alloc_grid = nullptr ; 
        } else {
            delete j_rgrid_f ; 
        }
        j_rgrid_f = nullptr ; 
    }


    List rgrid_b_final(len_rgrid_x) ; 
    for (int i = 0 ; i < len_rgrid_x  ; i++){
        NumericVector mat_b(pad_mat_rgrid_b[i].begin() , pad_mat_rgrid_b[i].end()) ; 
        mat_b.attr("dim") = IntegerVector(list_padded_dims[i].begin() , list_padded_dims[i].end()) ; 
        rgrid_b_final[i] = mat_b ; 
    }
    NumericVector simulation_final(Sim_prediction.begin() , Sim_prediction.end()) ; 
    List rgrid_x_final = wrap(pad_rgrid_x) ; 
    //pas besoin de wrap rgrid_f , elle est toujours sous le format Rcpp 
    List output_final = List::create(rgrid_x_final,rgrid_b_final,mat_rgrid_f,simulation_final) ; 

  
    return (output_final) ; 
}

