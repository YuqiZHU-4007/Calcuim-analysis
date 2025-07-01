


zi = 24;
ti = 100;

imio = ImageIO('Z:\20150609\rec7161475_sig\registered');
frame1 = readframe(imio, ti, zi);

imio = ImageIO('Z:\20150609\rec7161475_sig\dfdf_show_all');
frame2 = readframe(imio, ti, zi);

imio = ImageIO('Z:\20150609\rec7161475.dcimg', 45, 45);
frame3 = readframe(imio, ti, zi);

imio = ImageIO('Z:\20150609\rec7161475', 45, 45);
frame4 = readframe(imio, ti, zi);

sum(frame1(:)==frame3(:)) == numel(frame1)
%sum(frame2(:)==frame3(:)) == numel(frame1)
sum(frame3(:)==frame4(:)) == numel(frame1)

