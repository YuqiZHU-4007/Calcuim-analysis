dirf = getDirectory("Choose the Input Directory:functionnal");
savep = dirf;//substring(dirf,0,lengthOf(dira)-9);

list = getFileList(dirf); 
for (i=0; i<list.length; i++) {   
filename = dirf+list[i];
aname=list[i];
if (endsWith(filename, ".avi")) { 
	open(filename);fname=getTitle();
	fname=getTitle();
	run("Z Project...", "stop=43 projection=[Max Intensity]");
	aname=getTitle();
	saveAs("Tiff", savep+aname);
	aname=getTitle();
	//if (i==1) {nameb=aname };
	//if (i==list.length) {namea=aname };
	close(aname);
	close(fname);
	}
}
print("down");
nameb= dirf+"MAX_stim avr corr mapping37-42.tif";//'MAX_stim avr corr mapping37-42.tif';'MAX_stim avr corr mapping37-42_thr_0$65.tif'
namea= dirf+"MAX_stim avr corr mapping1-6.tif";
open(nameb);
open(namea);
//run("Merge Channels...", "c1=["+namea+"] c2=["+nameb+"]");
run("Merge Channels...", "c1=[MAX_stim avr corr mapping37-42.tif] c2=[MAX_stim avr corr mapping1-6.tif]");
saveAs("Tiff", savep+"merge_before_after.tif");
name=getTitle();
close(name);

