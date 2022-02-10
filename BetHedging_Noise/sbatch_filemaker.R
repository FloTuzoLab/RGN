param=read.table("/beegfs/home/flabourel/BetHedging/bin/IO_files/param.txt", h=TRUE)

lines=param$line

lines2=c()

for(i in lines){
    system("mkdir /beegfs/home/flabourel/BetHedging/bin/o")
    system("mkdir /beegfs/home/flabourel/BetHedging/bin/e")
    write(file=paste("simul", i, ".slurm", sep=""), c("#!/bin/csh", "#SBATCH -m n", "#SBATCH -p normal", "#SBATCH -N 1","#SBATCH --ntasks-per-node=1","#SBATCH --mem=500MB", "#SBATCH --time=24:00:00", "#SBATCH -o beegfs/home/flabourel/BetHedging/o", "#SBATCH -e /beegfs/home/flabourel/BetHedging/e", "cd /beegfs/home/flabourel/BetHedging/bin", paste("./simul ", i, sep="")))
    lines2=c(lines2, i)

}
write(file="run.sh", paste("sbatch simul", lines2, ".slurm", sep=""))
