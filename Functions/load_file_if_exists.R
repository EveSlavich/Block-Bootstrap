load_file_if_exists = function(path){
 if(file.exists(path)){
  print("file exists")
 
   load(path, envir=.GlobalEnv)
 }
  
}

