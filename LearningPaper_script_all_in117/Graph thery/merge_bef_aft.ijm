dirf = getDirectory("Choose the Input Directory:functionnal");
savep = dirf;//substring(dirf,0,lengthOf(dira)-9);

filename=dirf+"stim avr corr mapping1-6_auto_thr_aft.gif";
open(filename);
run("Z Project...", "stop=24 projection=[Max Intensity]");

filename=dirf+"stim avr corr mapping37-42_auto_thr_aft.gif";
open(filename);
run("Z Project...", "stop=24 projection=[Max Intensity]");

run("Merge Channels...", "c1=[MAX_stim avr corr mapping37-42_auto_thr_aft.tif] c2=[MAX_stim avr corr mapping1-6.tif]");
saveAs("Tiff", savep+"merge_before_after_auto_thr_aft.tif");