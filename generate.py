import os

list_mp = [500., 750., 1000., 1250., 1500.]
list_dl = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
condor_subs = []
model_dir = "/afs/ifh.de/group/atlas/users/pawlakd/modelsDM/DMRH"
mg_dir = "/afs/ifh.de/group/atlas/users/pawlakd/MG5_aMC_v3_5_5/bin/mg5_aMC"
working_dir = "/lustre/fs22/group/atlas/pawlakd/condor/"
log_dir = os.path.join(working_dir, "log")
os.makedirs(log_dir, exist_ok=True)
for mp in list_mp:
    for dl in list_dl:
        dirname = str(mp) + "_" + str(dl)
        out_path = os.path.join(log_dir, dirname+".out")
        os.system("touch %s"%out_path)
        log_path = os.path.join(log_dir, dirname+".log")
        os.system("touch %s"%log_path)
        err_path = os.path.join(log_dir, dirname+".err")
        os.system("touch %s"%err_path)
        dir_path = os.path.join(working_dir, dirname)
        os.makedirs(dir_path, exist_ok=True)
        script_name = os.path.join(dir_path, dirname + ".sh")
        submit_name = os.path.join(dir_path, dirname + "_sub")
        mg_com = os.path.join(dir_path, dirname + "_mg")
        condor_subs.append(submit_name)
        with open(script_name, "w") as script:
            script.write("#!/bin/bash\n")
            script.write("source ~/.bashrc\n")
            script.write(f"{mg_dir} {mg_com}\n")
        os.chmod(script_name, 0o777)  # Make the script executable
        with open(submit_name, "w") as submit:
            submit.write("universe = vanilla\n")
            submit.write("executable = %s\n"%script_name)
            submit.write('arguments  = \$(ClusterId) \$(ProcId)\n')
            submit.write(f"output = {working_dir}log/{dirname}.out\n")
            submit.write(f"error = {working_dir}log/{dirname}.err\n")
            submit.write(f"log = {working_dir}log/{dirname}.log\n")
            submit.write("request_memory=1024*20\n")
            submit.write("RequestCpus=4\n")
            submit.write('+JobFlavour = "tomorrow"\n')
            submit.write('queue 1')
        with open(mg_com, "w") as mg:
            mg.write(f"import model {model_dir}\n")
            mg.write("define chita = chit chit~\n")
            mg.write("define tta = t t~\n")
            mg.write("define phia = phi phi~\n")
            mg.write("define vl = ve vm vt\n")
            mg.write("define vl~ = ve~ vm~ vt~\n")
            mg.write("define l+ = e+ mu+ ta+\n")
            mg.write("define l- = e- mu- ta-\n")
            mg.write("define p = p b b~\n")
            mg.write("generate p p > phia tta chita\n")
            mg.write(f"output {dir_path}\n")
            mg.write("launch\n")
            mg.write("shower = Pythia8\n")
            mg.write("detector = Delphes\n")
            mg.write("madspin=on\n")
            mg.write("done\n")
            mg.write("/lustre/fs22/group/atlas/pawlakd/DMFV_RH_SFF_2med_nodecay/Cards/delphes_card_ATLAS.dat\n")
            mg.write("/lustre/fs22/group/atlas/pawlakd/DMFV_RH_SFF_2med_nodecay/Cards/pythia8_card.dat\n")
            mg.write("/lustre/fs22/group/atlas/pawlakd/DMFV_RH_SFF_2med_nodecay/Cards/run_card.dat\n")
            mg.write("/lustre/fs22/group/atlas/pawlakd/DMFV_RH_SFF_2med_nodecay/Cards/param_card.dat\n")
            mg.write("/afs/ifh.de/group/atlas/users/pawlakd/MG5_aMC_v3_5_5/bin/test2/Cards/madspin_card.dat\n")
            mg.write(f"set lam1 {dl}\n")
            mg.write(f"set lam2 {dl}\n")
            dl3 = dl + 1.0
            mg.write(f"set lam3 {dl3}\n")
            mg.write(f"set mphi {mp}\n")
            mg.write("set Wphi Auto\n")
            mg.write("done\n")

# Go back to the working directory
os.chdir(working_dir)

# Create the final submission script
script_mom = os.path.join(working_dir, "over.sh")
with open(script_mom, "w") as script:
    script.write("#!/bin/bash\n")
    for i in condor_subs:
        script.write(f"condor_submit {i}\n")

os.chmod(script_mom, 0o777)
