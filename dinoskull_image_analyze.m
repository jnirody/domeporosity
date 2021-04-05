 function dinoskull_image_analyze(directory,prefix,padding,type,start_slice,last_slice,voxel_size,normalize,calc_dist,background_black)

%% Function parameters (in order)
% directory, e.g., 'Data/AMNH5450_transverse'
% prefix, e.g., 'AMNH5450_hor_reslice'
% padding (how many zeros for single-digit numbers), e.g., '000'
% type ('dcm', 'jpg', or 'tif' (be careful about case!)), e.g., 'jpg'
% start_slice, e.g., 1
% last_slice, e.g., 100
% voxel_size (enter the voxel size in mm; enter )
% normalize (do you want to normalize slice numbers? (yes = 1, no = 0))
% calc_dist (do you want to calculate pore size stats? (yes = 1, no = 0))
% background_black (is the background black (yes = 1, no = 0))



%% read in dicom files and gather various information

if exist('background_black') == 0
    background_black = 1; 
end

if exist('calc_dist') == 0
    calc_dist = 0;
end

if exist('normalize') == 0
    normalize = 0;
end

add_on = 1-start_slice;


if strcmpi(type,'dcm') == 1
	dd = dicominfo(strcat(directory, prefix,padding(1:end-numel(num2str(start_slice))+1),num2str(start_slice)));
else
	dd = imfinfo(strcat(directory,prefix,padding(1:end-numel(num2str(start_slice))+1),num2str(start_slice),'.',type));
end
dimensions = [dd.Height dd.Width];	

%X_bw = zeros(dimensions(1),dimensions(2),last_slice-start_slice+1);
x_bw = zeros(dimensions(1),dimensions(2));

%X_filled = zeros(dimensions(1),dimensions(2),last_slice-start_slice+1);
x_fill = zeros(dimensions(1),dimensions(2));

total_area = zeros(last_slice-start_slice+1,1);
vasc_area = zeros(last_slice-start_slice+1,1);
if calc_dist == 1
    avg_area = zeros(last_slice-start_slice+1,1);
    std_area = zeros(last_slice-start_slice+1,1);
end

for slice = start_slice:last_slice
	take_off_pad = numel(num2str(slice)) - 1;
	filename = strcat(prefix,padding(1:end-take_off_pad),num2str(slice),'.',type);
	if strcmpi(type,'dcm') == 1
		image = dicomread(filename);
	else
		try
			image = imread(filename,type);
		catch 
			disp(strcat('Slice ',num2str(slice),' is missing.'));
            continue
		end
    end
	
	% binarize slice
	if background_black == 1
		level = graythresh(image);
        if slice < 200
            level = 0.25;
        end
		%X_bw(:,:,slice+add_on) = imcomplement(im2bw(image,level));
		x_bw = imcomplement(im2bw(image,level));
	else
		%X_bw(:,:,slice+add_on) = imcomplement(im2bw(image));
		x_bw = imcomplement(im2bw(image));
	end
	
 	%X_bw(:,:,slice+add_on) = bwdist(X_bw(:,:,slice+add_on)) <= 1;
	x_bw = bwdist(x_bw) <= 1;

	%X_bw(:,:,slice+add_on) = imcomplement(X_bw(:,:,slice+add_on));
	x_bw = imcomplement(x_bw);

 	%X_filled(:,:,slice+add_on) = imfill(X_bw(:,:,slice+add_on));
	x_fill(:,:) = imfill(x_bw,'holes');
	
 	%s1 = regionprops(X_filled(:,:,slice+add_on),'Area');
	s1 = regionprops(x_fill,'Area');
    
    if isempty(s1) == 1
        continue
    end

 	%s2 = regionprops(imcomplement(X_bw(:,:,slice+add_on) | ~X_filled(:,:,slice+add_on)),'Area'); 	
 	s2 = regionprops(~x_bw & x_fill,'Area');
	
    if slice == 200
        imshow(x_bw)
    end
	total_area(slice+add_on) = s1.Area*(voxel_size^2);
	vasc_area(slice+add_on) = sum(cat(1,s2.Area))*(voxel_size^2);
    disp(strcat('Slice ', num2str(slice)))
    
    if calc_dist == 1
        CC = bwconncomp(~x_bw & x_fill);
        CC_areas = regionprops(CC, 'Area');
        areas = zeros(length(CC_areas),1);
        for count = 1: length(CC_areas)
            areas(count) = CC_areas(count).Area*(voxel_size^2);
        end
        
        avg_area(slice+add_on) = mean(areas);
        std_area(slice+add_on) = std(areas);

    end
end
percent_vasc = vasc_area./total_area
avg_percent_vasc = mean(percent_vasc)

%temp = directory;
%outfile = directory;
%while temp
%	[outfile,temp] = strtok(temp,'/');
%end

outfile = prefix;

% Percent vascularity figure
if normalize == 1
    xaxis = linspace(0,1,last_slice-start_slice+1);
else
    xaxis = start_slice:last_slice;
end

figure; plot(xaxis,percent_vasc);
ylabel('Percent vascularity');
saveas(gcf,strcat(outfile,'_pervasc'),'epsc');

if calc_dist == 1
    % Average pore area figure
    figure; plot(xaxis,avg_area); 
    ylabel('Average pore area (mm^{2})');
    saveas(gcf,strcat(outfile,'_avgporearea'),'epsc');
    
    % Avg + Std pore area figure
    figure; shadedErrorBar(xaxis,avg_area,std_area); 
    ylabel('Average (Standard Deviation) pore area (mm^{2})');
end


if calc_dist == 1
    save(strcat(outfile,'.mat'),'total_area','vasc_area','percent_vasc','avg_percent_vasc','avg_area','std_area');
else
    save(strcat(outfile,'.mat'),'total_area','vasc_area','percent_vasc','avg_percent_vasc');
end

end

