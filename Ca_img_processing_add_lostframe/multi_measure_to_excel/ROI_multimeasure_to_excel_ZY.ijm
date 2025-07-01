dir = getDirectory("Choose the Input Directory");
dirresult=getDirectory("Choose the Output Directory");
dirroi=getDirectory("Choose the ROI Directory");
setBatchMode(true); 
list = getFileList(dir);
for (i=0; i<list.length; i++) {    
showProgress(i+1, list.length);
filename = dir+list[i];
roiname=list[i];
roiname =substring(list[i], 0, lengthOf(list[i])-1);
open(filename, "virtual");
Name=getInfo("image.filename");
roiManager("Open", dirroi+roiname+".zip");
roiManager("Multi Measure");
saveAs("Results", dirresult +roiname+".xls");
run("close table jru v1", "table=Results");
roiManager("Delete");
run("Close");
}