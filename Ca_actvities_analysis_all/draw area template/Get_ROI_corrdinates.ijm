dir = getDirectory("Choose the Input Directory");
list = getFileList(dir); 
for (i=0; i<list.length; i++) {    

filename = dir+list[i];
if (endsWith(filename, ".zip")) {  
	if (1) {
	//
	roiManager("Open", filename);
	print("roiname roi_point_No X Y Z");
	//get roi total number 
	var roinumber=roiManager("count");
	//get roi x/y coordinate
	for (roii = 0; roii < roinumber; roii++) {
		//select roi, and get the roi name
		roiManager("select",roii);
		roiname=Roi.getName;
		//print(roiname);
		//get coordinate
		Roi.getCoordinates(x, y); //edge point
		//get Contained Points
		//Roi.getContainedPoints(x, y); // area point
		for (roipi = 0; roipi < x.length; roipi++) {
			//output the x,y,z
			print(roiname,roipi,x[roipi], y[roipi],getSliceNumber());
		}				
	}
}
roiManager("deselect")
roiManager("Delete");
selectWindow("Log");
saveAs("Text", filename);
close("Log");
//Table.deleteRows(0, 0, "Log");
}
}
