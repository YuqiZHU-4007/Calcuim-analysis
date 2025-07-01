
dir = getDirectory("Choose the Input Directory");
dirresult=getDirectory("Choose the Output Directory");
dirroi=getDirectory("Choose the ROI Directory");
setBatchMode(true); 
list1 = getFileList(dir);
for (i=0; i<list1.length; i++) {    
showProgress(i+1, list1.length);
roiname=list1[i];roiname =substring(list1[i], 0, lengthOf(list1[i])-1);
roiManager("Open", dirroi+roiname+".zip");
filename = dir+list1[i];list = getFileList(filename);
for (ii=0; ii<list.length; ii++) { 
	ffilename=filename+list[ii]; 
	open(ffilename);  
	outputname=dirresult +list1[i]+substring(list[ii], 0, lengthOf(list[ii])-4)+".xls";
roiManager("Multi Measure");
saveAs("Results", outputname);
run("close table jru v1", "table=Results");
run("Close");

}
roiManager("Delete");
}
