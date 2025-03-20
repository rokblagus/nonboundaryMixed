

####scripts for lmm

  
zz=0
for (N in c(25,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.01,0.1)){
      for (mult in c(0.5,1.5)){
        for (rh in c(0.5,0.8)){
          mod="lmm"
    zz=zz+1        
      cat(
        "set.seed(",zz,")\n",
        "num_cluster =",N,"\n",
        "num_subj=",n,"\n",
        "sd_int=sqrt(",vr_int,")\n",
        "multiplier_int=",mult,"\n",
        "rho=",rh,"\n",
        "model=\"",mod,"\"\n",
        "source(\"../source_fun.R\")"
        ,file=paste0("lmm/script_",zz,".R"),sep="")
                  
            
          }
}}}}



####binom

zz=0
for (N in c(25,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.01,0.1)){
      for (mult in c(0.5,1.5)){
        for (rh in c(0.5,0.8)){
          mod="bin"
          zz=zz+1        
          cat(
            "set.seed(",zz,")\n",
            "num_cluster =",N,"\n",
            "num_subj=",n,"\n",
            "sd_int=sqrt(",vr_int,")\n",
            "multiplier_int=",mult,"\n",
            "rho=",rh,"\n",
            "model=\"",mod,"\"\n",
            "source(\"../source_fun.R\")"
            ,file=paste0("bin/script_",zz,".R"),sep="")
          
          
        }
      }}}}




####pois

zz=0
for (N in c(25,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.01,0.1)){
      for (mult in c(0.5,1.5)){
        for (rh in c(0.5,0.8)){
          mod="pois"
          zz=zz+1        
          cat(
            "set.seed(",zz,")\n",
            "num_cluster =",N,"\n",
            "num_subj=",n,"\n",
            "sd_int=sqrt(",vr_int,")\n",
            "multiplier_int=",mult,"\n",
            "rho=",rh,"\n",
            "model=\"",mod,"\"\n",
            "source(\"../source_fun.R\")"
            ,file=paste0("pois/script_",zz,".R"),sep="")
          
          
        }
      }}}}

