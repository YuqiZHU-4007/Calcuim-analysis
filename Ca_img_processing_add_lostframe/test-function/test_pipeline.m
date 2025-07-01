%test 配准程序（用tuborreg）
savepath=checkpath([opt.savepath '/test for 20190903/']);
imio = dcimgreader_zyq(opt.imgpath, opt.T, [],opt.process);  
slice = double(imio.readframe(833, 1));imwrite(uint16(slice), [opt.savepath '/' 'slice' '.tif']);%figure,imshow(uint16(slice),'InitialMagnification','fit');
refr_frame = env.vol(:,:,1); 

%算volmask和mask
[slice_volmask,slice_mask]=getmask(slice,0,opt,env);
imwrite(uint16(slice_mask), [savepath '/' 'slice_mask' '.tif']);
imwrite(uint16(slice_volmask), [savepath '/' 'slice_volmask' '.tif']);

[ref_volmask,ref_mask]=getmask(refr_frame,0,opt,env);
imwrite(uint16(ref_mask), [savepath '/' 'reg_mask' '.tif']);
imwrite(uint16(ref_volmask), [savepath '/' 'reg_volmask' '.tif']);

%去背景
slice_debackground=slice;slice_debackground(~slice_mask)=false;imwrite(uint16(slice_debackground), [savepath '/' 'slice_debackground' '.tif']);
ref_debackground=refr_frame;ref_debackground(~ref_mask)=false;
%配准
%配准去背景(包括模板) reg_frame1/2/3
pre_refr = prepare_reference_image_for_registration(ref_debackground, opt.warpopt);
[reg_frame, reg_params] = frame_registration_prepared_reference(pre_refr, slice_debackground, opt.warpopt);
imwrite(uint16(reg_frame), [savepath '/' 'reg_frame_m1' '.tif']);%figure,imshow(uint16(slice),'InitialMagnification','fit');
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
%去背景(不模板)
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
[reg_frame, reg_params] = frame_registration_prepared_reference(pre_refr, slice_debackground, opt.warpopt);
imwrite(uint16(reg_frame), [savepath '/' 'reg_frame_m2' '.tif']);%figure,imshow(uint16(slice),'InitialMagnification','fit');
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
%不去背景
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
[reg_frame, reg_params] = frame_registration_prepared_reference(pre_refr, slice, opt.warpopt);
imwrite(uint16(reg_frame), [savepath '/' 'reg_frame_m3' '.tif']);%figure,imshow(uint16(slice),'InitialMagnification','fit');

%去背景
slice_debackground=slice;slice_debackground(~slice_volmask)=false;
ref_debackground=refr_frame;ref_debackground(~ref_volmask)=false;
%配准
%配准去背景(包括模板) reg_frame1/2/3
pre_refr = prepare_reference_image_for_registration(ref_debackground, opt.warpopt);
[reg_frame, reg_params] = frame_registration_prepared_reference(pre_refr, slice_debackground, opt.warpopt);
imwrite(uint16(reg_frame), [savepath '/' 'reg_frame_vm1' '.tif']);%figure,imshow(uint16(slice),'InitialMagnification','fit');
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
%去背景(不模板)
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
[reg_frame, reg_params] = frame_registration_prepared_reference(pre_refr, slice_debackground, opt.warpopt);
imwrite(uint16(reg_frame), [savepath '/' 'reg_frame_vm2' '.tif']);%figure,imshow(uint16(slice),'InitialMagnification','fit');
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
%不去背景
pre_refr = prepare_reference_image_for_registration(refr_frame, opt.warpopt);
[reg_frame, reg_params] = frame_registration_prepared_reference(pre_refr, slice, opt.warpopt);
imwrite(uint16(reg_frame), [savepath '/' 'reg_frame_vm3' '.tif']);%figure,imshow(uint16(slice),'InitialMagnification','fit');