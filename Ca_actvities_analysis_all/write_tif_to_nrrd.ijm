dir = getDirectory("Choose the Input Directory");
list = getFileList(dir); 
for (i=0; i<list.length; i++) {    

filename = dir+list[i];
Cname=list[i];
if (endsWith(filename, ".tiff")) {  
open(filename);
Name=substring(filename,0,lengthOf(filename)-4);
run("Make Substack...", "  slices=1-24");
filename=getTitle();
run("Properties...", "channels=1 slices=24 frames=1 unit=micro pixel_width=0.66 pixel_height=0.66 voxel_depth=8");
run("Nrrd ... ", "nrrd=["+Name+".nrrd]");

run("Close All");
 } 
}