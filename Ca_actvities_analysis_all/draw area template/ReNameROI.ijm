
//rename roi, for import roi info
if (0) {
	var rois=roiManager("count");
	for (i = 0; i < rois; i++) {
		roiManager("select", i);
		name=Roi.getName;
		roiManager("rename", name);
	}
}

//name roi
if (0) {
	
	//you can change the str "cell#" to any you like str
	
	//get roi total number 
	var roinumber=roiManager("count");
	print(roinumber);
	var cellnumber=0;
	var igvalue=0;
	var oldcellnumber=0;
	for (i = 0; i < roinumber; i++) {
		roiManager("select", i);
		roiname=Roi.getName;
	
		//check the ignore label	
		if (roiname=="ignore") {igvalue=1;continue;};
		//chcek the new cell roi after the ignore
		if (lengthOf(roiname)<5 && roiname!="ignore") {igvalue=0;};
	
		//rename ignore roi
		if (igvalue==1 && lengthOf(roiname)==14) {roiname="ignore";};
	
		//skip had rename roi
		namelist=split(roiname,"#");
		if (namelist[0]=="cell") {
			cellnumber=parseInt(namelist[1]);
			if(cellnumber>oldcellnumber){oldcellnumber=cellnumber;};		
			continue;
			};
				
		//rename roi	
		if (lengthOf(roiname)==14 && igvalue==0) {			
			roiname="cell#"+cellnumber;
			}else if (igvalue==0 && roiname==0) {
				if(oldcellnumber!=0){cellnumber=oldcellnumber;
				};			
				cellnumber++;
				oldcellnumber=cellnumber;
				roiname="cell#"+cellnumber;
				}else if(igvalue==0 && roiname!=0){
				oldcellnumber=cellnumber;
				cellnumber=parseInt(roiname);
				roiname="cell#"+cellnumber;			
				};	
		
		//print(cellnumber);
		roiManager("rename",roiname);
		print(Roi.getName);			
		};	
}
	
//save image
if (0) {
//save image
//run("From ROI Manager");
//var figname=split(getTitle(), ".");
//saveAs("Tiff", getDirectory("Choose a Directory")+"/"+figname[0]);		
}


//get ROI range Coordinate
if (0) {
	//
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
		//Roi.getCoordinates(x, y);
		//get Contained Points
		Roi.getContainedPoints(x, y);
		for (roipi = 0; roipi < x.length; roipi++) {
			//output the x,y,z
			print(roiname,roipi,x[roipi], y[roipi],getSliceNumber());
		}				
	}
}

