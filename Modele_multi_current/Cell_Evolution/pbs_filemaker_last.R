param=read.table("/panhome/flabourel/Cell_Evo_last/bin/IO_files/param.txt", h=TRUE)

lines=param$line

lines2=c()

for(i in lines){
    system("mkdir /panhome/flabourel/Cell_Evo_last/bin/o")
    system("mkdir /panhome/flabourel/Cell_Evo_last/bin/e")
    write(file=paste("simullast", i, ".pbs", sep=""), c("#!/bin/csh", "#PBS -m n", "#PBS -q q1week", "#PBS -l nodes=1:ppn=1,mem=1gb", "#PBS -o /panhome/flabourel/Cell_Evo_last/o", "#PBS -e /panhome/flabourel/Cell_Evo_last/e", "cd /panhome/flabourel/Cell_Evo_last/bin", paste("./simullast ", i, sep="")))
    lines2=c(lines2, i)

}
write(file="run_last.sh", paste("qsub simullast", lines2, ".pbs", sep=""))
