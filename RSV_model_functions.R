# RSV_model_functions
fun_sub2ind=function(i_inf,j_age,varname,varname_list,n_var,n_age,n_inf){
  varnum=which(varname_list %in% varname); k=(j_age-1)*n_var*n_inf + (varnum-1)*n_inf + i_inf; k }

fun_sirs_varnames=function(varname_list,n_age,n_inf){
  array(sapply(1:n_age, function(x_age) {sapply(varname_list, function(x) {paste0(paste0(x,'_',1:n_inf),'_',x_age)})} ))
}