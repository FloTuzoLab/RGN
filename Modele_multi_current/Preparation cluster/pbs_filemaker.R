param=read.table("mkdir /panhome/flabourel/Cell_Evolution/bin/IO_files/param.txt", h=TRUE)

lines=param$line

lines2=c()

for(i in lines){
    system("mkdir /panhome/flabourel/Cell_Evolution/bin/o")
    system("mkdir /panhome/flabourel/Cell_Evolution/bin/e")
    write(file=paste("simul", i, ".pbs", sep=""), c("#!/bin/csh", "#PBS -m n", "#PBS -q q1week", "#PBS -l nodes=1:ppn=1,mem=1gb", "#PBS -o /panhome/flabourel/Cell_Evolution/o", "#PBS -e /panhome/flabourel/Cell_Evolution/e", "cd /panhome/flabourel/Cell_Evolution/bin", paste("./simul ", i, sep="")))
    lines2=c(lines2, i)

}
write(file="run.sh", paste("qsub simul", lines2, ".pbs", sep=""))
