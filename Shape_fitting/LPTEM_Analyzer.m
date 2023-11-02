%% Liquid Phase Electron Microscopy Video Analyzer
% Reads in-situ video files with the .dm4 format
% Stores relevant information in the video_info structure
% Written by Ivan A. Moreno-Hernandez, starting on 02/24/2022

%Clear all variables
clear all
clc

% select folder that contains several video folders
number_folders = 9 % Select number of folders to process
for i = 1:number_folders
i %indicates which folder to select
folders{i}= uigetdir('','*Select Folder with video sub-folders*'); %Read directory
end
%% Create video_info structure with information about videos
count = 1;
for i = 1:number_folders
clear video_info_temp
%Obtain information about videos in temporary variable
video_info_temp = prepare_directory(folders{i});
video_info_temp= last_frame_info(video_info_temp);
video_info_temp = find_true_fps(video_info_temp);
%Store all information from all videos in video_info structure
for  j = 1:size(video_info_temp,2)
video_info(count) = video_info_temp(j); %Construct 
count = count +1;
end

end
clearvars -except video_info %Clear all variables except video_info

save RuO2_manuscript_video_info video_info
clear all
%% Load previous video_info and find start of frames
clear all
clc
load  RuO2_manuscript_video_info video_info
% Find start of videos in chunks
video_info = find_start(video_info,1:102); %This will take about an hour

%% Define the time domain for all videos
clear all
load  RuO2_manuscript_video_info video_info
video_info = lc_time(video_info,1:size(video_info,2));

%% Fix intensity and resolution scales for videos that were not read properly
for i = 1:size(video_info,2)
    if video_info(i).scale_nm_pix == 1
    video_info(i).scale_nm_pix = 0.21252881; %For 2k videos collected in this study
    end

    if video_info(i).intscale == 1
    video_info(i).intscale = 0.018656718; %For Duke's OneView
    end

end

save  RuO2_manuscript_video_info video_info
clear all
%% %% Use this section to view videos and make a comment about viability for analysis
clear all
load RuO2_manuscript_video_info video_info
%%
for i = 100
v = i
sp = 64; %Speed factor, greater than 1
video_info = playback_video(video_info,v,sp)
video_info(v).analysis = questdlg('Is the video suitable for further analysis?');
end
%% %2023/09/11: start here
clear all

load  RuO2_manuscript_video_info video_info
%% Use this section to define lines of interest for each nanocrystal frame

v = 102;
if isequal(video_info(v).analysis,'Yes')
video_info = loi_gui(video_info,v,15)
end


%% Pull frames from video
%Show first and last frames of video
%Define image parameters
I_level = 1; %Brightness level, 0-1
aspect = 2; %Aspect ratio of image
im_int = 5; %Frames used for image integration (integer, frames +- the frames listed)
end_factor = .85;
% Change number of frames to define resolution of video)
n_frames = 20; %Number of frames to analyze
end_frame = round(end_factor*size(video_info(v).frame_list,1));
frames = round(linspace(video_info(v).start_frame,end_frame,n_frames));

%Pull line infro from_video_info and define aspect ratio
line_info_all = video_info(v).lines;
line_length = max(round(((line_info_all(:,1)-line_info_all(:,3)).^2+(line_info_all(:,2)-line_info_all(:,4)).^2).^(1/2)))
h = round(line_length*1.3); %Height of box in pixels
w = round(h/aspect);  %Width of box in pixels

I_start = crop_box_line2(video_info,v,frames(1),h,w,I_level,im_int);
I_end = crop_box_line2(video_info,v,frames(end),h,w,I_level,im_int);
imshowpair(I_start,I_end,'montage')


% Save parameters
%save all the info to video_info
video_info(v).h = h;
video_info(v).w = w;
video_info(v).I_level = I_level;
video_info(v).im_int = im_int;
video_info(v).frames = frames;



%% Segment the videos by hand
for i = 1:1:size(video_info(v).frames,2)
complete = i/size(video_info(v).frames,2)*100
    
     I = crop_box_line2(video_info,v,video_info(v).frames(i),video_info(v).h,video_info(v).w,.5,1); % Get image
     I = imadjust(I);
     imshow(I)
    roi = drawpolygon
    Ibw = createMask(roi);
    imshow(Ibw)
    video_info(v).Ibw{video_info(v).frames(i)} = Ibw;
    clf

 end

save  RuO2_manuscript_video_info video_info


%% Prepare initial guess
    L_guess(1) = 450;
    L_guess(7) = L_guess(1)*6.6; %Scale length
    L_guess(5) = 170;
    L_guess(6) = 200;
    L_guess(3) = L_guess(7)*0.84; %Around 0.75 default
    L_guess(4) = L_guess(7)*0.71; %Around 0.75 default
    angles = [45+00, 90, 90-0];
    translation_shift = [-0,-0]; %Change translation as necessary

 %Keep this one the same
  L_guess(2) = L_guess(1);
 pixel_scale_A = double(video_info(v).scale_nm_pix*10); %A/pixel
 translation = [video_info(v).h/2,video_info(v).w/2]*pixel_scale_A;
translation = translation + translation_shift;
 
 [~,shp2d] = nanocrystal_RuO2_shape(L_guess,angles,translation);
     Ibw =  video_info(v).Ibw{video_info(v).frames(1)};

    Ibw_test = ibw_shape(Ibw,pixel_scale_A,shp2d);
    imshowpair(Ibw,Ibw_test)
video_info(v).L_guess_initial = L_guess;


%% Save guess %End here
save  RuO2_manuscript_video_info video_info

%% load data
clear all
load  RuO2_manuscript_video_info video_info

%% Fit frames

clc
for v = 3
    v
try
clear fit_all
try
fit_all = video_info(v).fit_all;
end

for i = 1:size(video_info(v).frames,2) %make par'
   
test = 0;
    try
 if fit_all(i,1) > 0 %Data already fit
     test = 1;
     
end
    end
 
 try

     if test ~= 1;
         [v,i]
Ibw =  video_info(v).Ibw{video_info(v).frames(i)};
imshow(Ibw)
%Calculate picture limits
lim_x = [min(find(sum(Ibw,2) > 0))-50,max(find(sum(Ibw,2) > 0))+50];
lim_y = [min(find(sum(Ibw,1) > 0))-50,max(find(sum(Ibw,1) > 0))+50];

%Limiters
lim_x = max(lim_x,1);
lim_y = max(lim_y,1);
lim_x = min(lim_x,size(Ibw,1));
lim_y = min(lim_y,size(Ibw,2));

%Make smaller image to decrease computational cost
Ibw2 = Ibw(lim_x(1):lim_x(2),lim_y(1):lim_y(2));

imshow(Ibw2)
L_guess = video_info(v).L_guess_initial*.1; %Use 0.5 multiplier to fix guesses
pixel_scale_A = double(video_info(v).scale_nm_pix*10); %A/pixel
translation = [size(Ibw2,1)/2,size(Ibw2,2)/2]*pixel_scale_A;
%translation = [video_info(v).h/2,video_info(v).w/2]*pixel_scale_A;
angles = [45, 90, 90];
fit_scale_down = 2^0;
tic
%[L_fit,angle_fit,translation_fit] = fit_RuO2_v3(L_guess,angles,translation,imresize(double(Ibw),1/fit_scale_down),pixel_scale_A*fit_scale_down)
[L_fit,angle_fit,translation_fit] = fit_RuO2_v3(L_guess,angles,translation,imresize(double(Ibw2),1/fit_scale_down),pixel_scale_A*fit_scale_down)

fit_all(i,:) = [L_fit,angle_fit,translation_fit];

i
     end
 end
video_info(v).fit_all = fit_all;

end

end
end

save  RuO2_manuscript_video_info video_info

%% View fits
clear all
load RuO2_manuscript_video_info video_info

%%
clc
for v = 102
try
for i = 1:size(video_info(v).frames,2) %make par for higher throughput
[v,i]
Ibw =  video_info(v).Ibw{video_info(v).frames(i)};
imshow(Ibw);
%Calculate picture limits
lim_x = [min(find(sum(Ibw,2) > 0))-10,max(find(sum(Ibw,2) > 0))+10];
lim_y = [min(find(sum(Ibw,1) > 0))-10,max(find(sum(Ibw,1) > 0))+10];

%Limiters
lim_x = max(lim_x,1);
lim_y = max(lim_y,1);
lim_x = min(lim_x,size(Ibw,1));
lim_y = min(lim_y,size(Ibw,2));

%Make smaller image to decrease computational cost
Ibw2 = Ibw(lim_x(1):lim_x(2),lim_y(1):lim_y(2));
fit_all=video_info(v).fit_all;
[~,shp2d] = nanocrystal_RuO2_shape(fit_all(i,1:7),fit_all(i,8:10),fit_all(i,11:12));
pixel_scale_A = double(video_info(v).scale_nm_pix*10); %A/pixel
Ibw_fit = ibw_shape(Ibw2,pixel_scale_A,shp2d);
clf
tiledlayout(1,2)
nexttile
imshowpair(Ibw2,Ibw_fit)
nexttile
hold on
plot(fit_all(:,1))
plot(i,fit_all(i,1),'o')

hold off
set(gcf,'Position',[1500 10 1000 1000])
drawnow
pause(.1)
end


end
end

%% Count good videos

good_total = 0;
total_videos = 0;
for i = 1:102
try
    if isequal(video_info(i).comment,'Good')
    good_total = good_total +1;
    video_info_good(good_total) = video_info(i);
    end
        if isequal(video_info(i).analysis,'Yes')
    total_videos = total_videos +1;
    end
end
end
save  RuO2_manuscript_video_info_good video_info_good

%% Functions for Script


% This part of the code is used to prepare the locations of all of the
% frames from the videos
function out = prepare_directory(x)
folders = folders_in(x);
for i = 1:size(folders,1)
    out(i).video_folder = folders(i);
    out(i).frame_list = find_frames(folders(i));
    file_dir = "";
    name = "";
    for j = 1:size(out(i).frame_list,1)
        file_dir = [file_dir; strcat(out(i).frame_list(j).folder,'\',out(i).frame_list(j).name)];
        name = [name; out(i).frame_list(j).name];
    end
    file_dir = [file_dir(2:end)]; % Remove initial empty string
    name = [name(2:end)]; % Remove initial empty string
    out(i).file_dir = file_dir;
    out(i).name = name;
    out(i).num_frames = size(file_dir,1);
    
end
end


% Find all folders within a folder
function out = folders_in(x)
d = dir(x);
% remove all files (isdir property is 0)
dfolders = d([d(:).isdir]==1);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
out = strings(length(dfolders),1);
for floop = 1:length(dfolders)
out(floop,1) = strcat(dfolders(floop).folder,'\',dfolders(floop).name);
end
end


% find all .dm3 frames in all folders and subfolders
function out = find_frames(x)
filelist = dir(fullfile(x, '**\*.dm4'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list
out = filelist;
end


% Determines information from last video frame
function x = last_frame_info(x)

for i = 1:size(x,2)
    i
    Imagestruc = dmread(x(i).file_dir(end)); %reads dm4 file
    x(i).x_pix = Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed0.Value; %extracts x-dimension of image
    x(i).y_pix =Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed1.Value; %extracts y-dimension of image
    x(i).scale_nm_pix =Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Dimension.Unnamed0.Scale.Value; %extracts scale of image (nm/pixel)
    x(i).intscale =Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Brightness.Scale.Value ; %extracts intensity scale for image (e-/value)
    x(i).fps = 1/Imagestruc.ImageList.Unnamed0.ImageTags.DataBar.ExposureTime0x28s0x29.Value; % FPS = 1/exposure
    x(i).final_dose =  x(i).fps*mean(Imagestruc.ImageList.Unnamed0.ImageData.Data.Value)*Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Brightness.Scale.Value/((100)*x(i).scale_nm_pix.^2); % image dose (e- A^-2 s^-1)
end
end



%Find the frames per second for the video
function x = find_true_fps(x)
for i = 1:size(x,2)
file_list =[];
file_list = vertcat(x(i).frame_list(:).folder);
count = 0;
for j = 1:size(file_list,2)
    if file_list(j,:) == file_list(1,:)
    count = count +1
    end
end
x(i).fps = count;
end
end

%Find beginning of video based on intensity
function x  = find_start(x,k)
parfor i = k %Can switch to par for for speed
    i
%Determine initial frame by comparing electrons detected to previous frames
count = size(x(i).file_dir,1);
trajectory_dose_frac =1;
%Fast backward scan
I_final_sum = sum(read_frame(x(i).file_dir(end)),'all');
while trajectory_dose_frac > 0.5; 
    Im = read_frame(x(i).file_dir(count)); % Read frame from count
    trajectory_dose_frac = sum(Im,'all')./I_final_sum;
 count = count-10%
 if count < 1
     count = 1;
     trajectory_dose_frac =0
 end
end

%Frame-by-frame forwards scan
while trajectory_dose_frac < 0.7;
if count < 1
     count = 1;
end
    Im = read_frame(x(i).file_dir(count)); % Read frame from count
    trajectory_dose_frac = sum(Im,'all')./I_final_sum;
 count = count+1 %
if count < 1
     count = 1;
end
end


x(i).start_frame = count;    
  
    
end
end



%reads dm4 file and outputs array in e-/(pixel)
function out = read_frame(x) %
    Imagestruc = dmread(x); 
    x_pix = Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed0.Value; %extracts x-dimension of image
    y_pix =Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed1.Value; %extracts y-dimension of image
    intscale =Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Brightness.Scale.Value ; %extracts intensity scale for image (e-/value)
    out = reshape(Imagestruc.ImageList.Unnamed0.ImageData.Data.Value,[x_pix,y_pix])'; %Normalize later
end



function x = lc_time(x,v)
for i = v
    
    x(i).time = (single([1:x(i).num_frames]') - single(x(i).start_frame))./single(x(i).fps);
end

end


%Watch video with imshow
function out = playback_video(x,v,sp)

out = x; %Make sure to keep data
frames_index = round(linspace(1,size(x(v).file_dir,1),round(size(x(v).file_dir,1)/sp)));
close all
f = figure(1)

for i = frames_index
I = read_frame(x(v).file_dir(i));
out(v).counts(i,:) = mean(I,"all");
I = I./max(I,[],'all');
I = imadjust(I);
title(['Frame: ', num2str(i)])
imshow(I);
f.Position = [1400 100 800 800];

end
end

%Line-of-interest graphical user interface
function x = loi_gui(x,i,n)

for j = i
framestart = x(j).start_frame;
frameend = x(j).num_frames;
frames = round(linspace(framestart,frameend,n));
for k = 1:n
I = read_frame(x(j).file_dir(frames(k)));
set(gcf,'name',x(j).name(frames(k)),'numbertitle','off')
imshow(imadjust(mat2gray(I)))
loi = drawline('InteractionsAllowed','all')
line_pos = customWait(loi)
line_info(k,1:4) = [loi.Position(1,:),loi.Position(2,:)]; %[x1 y1 x2 y2]?
line_tip(k,1:2) = [line_info(k,1),line_info(k,2)];
line_angle(k,1) = atan2(line_info(k,4)-line_info(k,2),line_info(k,3)-line_info(k,1))*(180/pi);

end

line_info_interp = interp1(frames,line_info,framestart:frameend);
line_tip_interp = interp1(frames,line_tip,framestart:frameend);
line_angle_interp = interp1(frames,line_angle,framestart:frameend)';

% Fill in for all frames prior to the start frame 
line_info_interp_all = uint32([repelem(line_info_interp(1,:),framestart-1,1); line_info_interp]);
line_info_interp_all = single(line_info_interp_all);

line_tip_interp_all = uint32([repelem(line_tip_interp(1,:),framestart-1,1); line_tip_interp]);

line_angle_interp_all = ([repelem(line_angle_interp(1,:),framestart-1,1); line_angle_interp]);
x(j).lines = line_info_interp_all;
x(j).lines_angle = line_angle_interp_all;
x(j).lines_tip = line_tip_interp_all;
end
close
end


% Wait function
function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;

end

function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end



function Iout = crop_box_line2(x,v,f,h,w,m,s)
% x = video_info
% v = video numbe 
% f = frame number
% h = height
% w = width
% m = mean value for image intensity
% s = number of frames +- the centered frame
I_all = [];
for j = [f-s:f+s]
    try
if isempty(I_all) == 1
I_all = read_frame(x(v).file_dir(j));
else
    I_all = I_all +read_frame(x(v).file_dir(j));
end

    end
end
I = I_all; %read_frame(x(v).file_dir(f));
I = I/mean(I,'all')*m;
[X,Y] = meshgrid(-(w-1)/2:(w-1)/2,-(h-1)/2:(h-1)/2);
XY = [X(:),Y(:)];
rot_mat = rotz(x(v).lines_angle(f));
XY_new = XY*rot_mat(1:2,1:2); %Rotate
line_info = x(v).lines(f,:);
line_center = [line_info(1)+line_info(3),line_info(2)+line_info(4)]/2;
XY_new = XY_new + [line_center(2),line_center(1)];
XY_new = round(XY_new);
size_I = size(I);
clear Intensity
for i= 1:size(XY_new,1)
    if XY_new(i,2) > size(I,2) | XY_new(i,2) < 1 | XY_new(i,1) > size(I,1) | XY_new(i,1) < 1 
    Intensity(i)  = m;
    else
    Intensity(i)  = I(XY_new(i,1),XY_new(i,2));
    end

end
Iout = flip(reshape(Intensity,size(X)),2);

end



function [shp3d,shp2d] = nanocrystal_RuO2_shape(L,angles,translation)
%RuO2 parameters
spacings = [4.4919, 4.4919,3.10660]; % In Angstroms %[4.4919, 4.4919, 3.1066];
unit_vectors = [1, 0, 0; 0, 1,0; 0,0,1];
b =  unit_vectors.*spacings';

%RuO2 Planes
L1 = L(1); %Side 1
L2 = L(2); %Side 2
L3 = L(3); %(111) tip 1 top
L4 = L(4); %(111) tip 2 top
L5 = L(5); %(111) tip 1 bottom
L6 = L(6); %(111) tip 2 bottom
L7 = L(7); %(001) tips
a = 1; %Integer values for for facets (aac)
c = 1; %Integer values for for facets (aac)
plane_info = [
    a, a, c, L3;
    -a, a, c, L3;
    a, -a, c, L4;
    -a, -a, c, L4;
    a, a, -c, L5;
    -a, a, -c, L5;
    a, -a, -c, L6;
    -a, -a, -c, L6;
    1, 1, 0, L1;
    -1, 1, 0, L2;
    -1, -1, 0, L1;
    1, -1, 0, L2;
    0, 0, 1, L7;
    0, 0, -1, 0];

[vertices_nanocrystal] = nanocrystal_shape2(plane_info,b);
vertices_nanocrystal = unique(vertices_nanocrystal,'rows'); %unique points
vertices_center = (max(vertices_nanocrystal)+min(vertices_nanocrystal))/2;
vertices_rotated = (vertices_nanocrystal-vertices_center)*rotz(angles(1))*roty(angles(3))*rotx(angles(2))+vertices_center;
vertices_rotated_translated = vertices_rotated + [translation,0]; %Translate particle
vertices_rotated_translated = double(vertices_rotated_translated); %Fix double error
shp3d = alphaShape(vertices_rotated_translated(:,1),vertices_rotated_translated(:,2),vertices_rotated_translated(:,3),Inf); %Output shape object
k = boundary(vertices_rotated_translated(:,1),vertices_rotated_translated(:,2),0); %indices for boundary
XY_boundary = [vertices_rotated_translated(k,1),vertices_rotated_translated(k,2)];
XY_boundary = unique(XY_boundary,'rows');
shp2d = alphaShape(XY_boundary(:,1),XY_boundary(:,2),Inf);
end




function Ibw_out = ibw_shape(I,pixel_scale_A,shp2d)
[X,Y] = meshgrid(1:size(I,1),1:size(I,2));  %Define Coordinates
X = X.*pixel_scale_A; %In angstroms
Y = Y.*pixel_scale_A; %In angstroms
I_inshape = inShape(shp2d,X(:),Y(:));
Ibw_out = reshape(I_inshape,size(X))';

end



function [L_fit,angle_fit,translation_fit] = fit_RuO2_v3(L,angles,translation,Icorr,pixel_scale_A)

L_fit = L;
angle_fit = angles;
translation_fit = translation;
for i = 1:1

[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,angle_fit] = dimension_fit_001_tip_global_rotation(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A)
%[L_fit,angle_fit] = dimension_fit_001_tip_global_rotation_fine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A)

[L_fit] = dimension_fit_reshape_111_tip(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_reshape_111_tip_part2(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);

[L_fit] = dimension_fit_110_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_111_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);

[L_fit] = dimension_fit_110_tip_global_refine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_111_tip_global_refine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_tip_global_refine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,translation_fit] = dimension_fit_001_pos(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,translation_fit] = dimension_fit_001_pos_v2(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_reshape_111_tip(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,angle_fit] = dimension_fit_001_tip_global_rotation(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A)
%[L_fit,angle_fit] = dimension_fit_001_tip_global_rotation_fine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A)
[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);

end

end



function [L_out] = dimension_fit_001_110(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);


[L110,L001] = meshgrid(-1:1,-1:1); %[110, 001]
step_size = 8;
fit_test = 0;
count_fits =  0;
while count_fits < 10; %10
while fit_test == 0
L_range = flip([L110(:),L001(:)]*step_size); %Flip fixes out of boundary search
for i = 1:9 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
L_loop(1) = L_loop(1)*(100+L_range(i,1))/100;
L_loop(2) = L_loop(2)*(100+L_range(i,1))/100;
L_loop(7) = L_loop(7)*(100+L_range(i,2))/100;
L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
L_trial(7) = L_trial(7)*(100+L_range(maxindx,2))/100;
L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;

end



%Computed similarity between images
function [sim_out] = im_similarity(Im_test,Im_actual)
% This function will calculate the similarity of an image, higher value is better
%Im_actual = imadjust(Im_actual);
%Im_actual = ordfilt2(Im_actual,1,ones(3,3));
Im_test  = imadjust(double(Im_test));
imshowpair(Im_test,Im_actual)
drawnow
try
    %im_diff = Im_test - Im_actual;
    
    im_diff = Im_test - Im_actual;
    t = im_diff >0; 
    im_diff(t) = 1.0*im_diff(t); %.71 to 1.0 for different error
    %(im_diff>0) = im_diff;
    %im_diff(im_diff<0) = im_diff*1;
    im_sum = sum(abs(im_diff),'all');
    im_sum;
   % sim_out = corr2(Im_test,Im_actual);
    sim_out = -im_sum; %1./(1+im_sum);
catch
    sim_out = -10^10;
end


end




%Reshapes the 111 tips
function [L_out] = dimension_fit_reshape_111_tip(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);


[Ltip1,Ltip2] = meshgrid(-1:1,-1:1); %[110, 001]
step_size = 100;
fit_test = 0;
count_fits =  0;
while count_fits < 30; %10
while fit_test == 0
L_range = [Ltip1(:),Ltip2(:)]*step_size;
for i = 1:9 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
%L_loop(3) = L_loop(3) + L_range(i,1);
L_loop(6) = L_loop(6) + L_range(i,1);
%L_loop(4) = L_loop(4) + L_range(i,2);
L_loop(5) = L_loop(5) + L_range(i,2);
%L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
%L_trial(3) = L_trial(3) + L_range(maxindx,1);
L_trial(6) = L_trial(6) + L_range(maxindx,1);
%L_trial(4) = L_trial(4) + L_range(maxindx,2);
L_trial(5) = L_trial(5) + L_range(maxindx,2);


%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
%L_trial(7) = L_trial(7)*(100+L_range(maxindx,2))/100;
%L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;

end

%
function [L_out] = dimension_fit_reshape_111_tip_part2(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);


[Ltip1,Ltip2] = meshgrid(-1:1,-1:1); %[110, 001]
step_size = 100;
fit_test = 0;
count_fits =  0;
while count_fits < 30; %10
while fit_test == 0
L_range = [Ltip1(:),Ltip2(:)]*step_size;
for i = 1:9 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
L_loop(3) = L_loop(3) - L_range(i,1);
L_loop(6) = L_loop(6) + L_range(i,1);
L_loop(4) = L_loop(4) - L_range(i,2);
L_loop(5) = L_loop(5) + L_range(i,2);
%L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
L_trial(3) = L_trial(3) - L_range(maxindx,1);
L_trial(6) = L_trial(6) + L_range(maxindx,1);
L_trial(4) = L_trial(4) - L_range(maxindx,2);
L_trial(5) = L_trial(5) + L_range(maxindx,2);


%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
%L_trial(7) = L_trial(7)*(100+L_range(maxindx,2))/100;
%L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;

end


%% Fit rotation and translation
function [angles_out,translation_out] = spatial_fit(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);


[X,Y,Z] = meshgrid(-1:1,-1:1,-1:1); %[X,Y,angle]

X = X;
Y = Y;
Z = Z/4; %Less deviation in angle
step_size = 16;
fit_test = 0;
count_fits =  0;
while count_fits < 6; %10
while fit_test == 0
XYZ_range = [X(:),Y(:),Z(:)]*step_size;
for i = 1:27 %parfor for speed
%Determine trial translation and angles
translation_loop = translation_trial +[XYZ_range(i,1),XYZ_range(i,2)];
angles_loop = angles_trial + [0,0,XYZ_range(i,3)];
%Compute shapes
try
[~,shp2d] = nanocrystal_RuO2_shape(L,angles_loop,translation_loop);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Ibw(end/2:end,:) =0; %
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
%drawnow
end
[~,maxindx] = max(corr_trial); %find best guess

translation_trial = translation_trial +[XYZ_range(maxindx,1),XYZ_range(maxindx,2)]; %Replaces old guess
angles_trial = angles_trial + [0,0,XYZ_range(maxindx,3)]; %Replaces old guess
fit_test = isequal(XYZ_range(maxindx,:),[0,0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
angles_out = angles_trial;
translation_out = translation_trial;


end



%Fits (111) tip ratios
function [L_out] = dimension_fit_111_tip_global_refine(L,angles,translation,I,pixel_scale_A)
L_trial = L;
L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_111_ratios = linspace(.98,1.02,51); %Ratios to use
for j = 3:6
    if j > 4
    Tip_111_ratios = linspace(.98,1.02,51); %Ratios to use

    end
for i = 1:size(Tip_111_ratios,2)
    
L_loop =  L_trial;
L_loop(j) =L_loop(j)*Tip_111_ratios(i); %Scale tip
    try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
    catch
corr_trial(i) = 0;
    end

%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(j) = L_trial(j)*Tip_111_ratios(maxindx);
end
L_out = L_trial;

end



%Fits (001) tip
function [L_out] = dimension_fit_001_tip_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_001_ratios = linspace(.5,1.5,141); %Ratios to use
for i = 1:size(Tip_001_ratios,2)
L_loop =  L_trial;
L_loop(7) =L_trial(7)*Tip_001_ratios(i); %Scale tip
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(7) = L_trial(7)*Tip_001_ratios(maxindx);

L_out = L_trial;

end


%Fits (001) tip
function [L_out] = dimension_fit_001_tip_global_refine(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_001_ratios = linspace(.98,1.02,51); %Ratios to use
for i = 1:size(Tip_001_ratios,2)
L_loop =  L_trial;
L_loop(7) =L_loop(7)*Tip_001_ratios(i); %Scale tip
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(7) = L_trial(7)*Tip_001_ratios(maxindx);

L_out = L_trial;

end






function [L_out] = dimension_fit_110_tip_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_110_ratios = linspace(.5,1.5,141); %Ratios to use
for i = 1:size(Tip_110_ratios,2)
L_loop =  L_trial;
L_loop(1) =L_trial(1)*Tip_110_ratios(i); %Scale sides
L_loop(2) =L_trial(2)*Tip_110_ratios(i); %Scale sides

try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(1) = L_trial(1)*Tip_110_ratios(maxindx);
L_trial(2) = L_trial(2)*Tip_110_ratios(maxindx);

L_out = L_trial;

end



%Fits (111) tip ratios
function [L_out] = dimension_fit_111_tip_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_111_ratios = linspace(.5,1.5,101); %Ratios to use
for j = 3:6
    if j > 4
    Tip_111_ratios = linspace(-.5,0.5,101); %Ratios to use

    end
for i = 1:size(Tip_111_ratios,2)
    
L_loop =  L_trial;
L_loop(j) =L_trial(7)*Tip_111_ratios(i); %Scale tip
    try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
    catch
corr_trial(i) = 0;
    end

%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(j) = L_trial(7)*Tip_111_ratios(maxindx);
end
L_out = L_trial;

end

%


function [L_out] = dimension_fit_110_tip_global_refine(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_110_ratios = linspace(.98,1.02,51); %Ratios to use
for i = 1:size(Tip_110_ratios,2)
L_loop =  L_trial;
L_loop(1) =L_loop(1)*Tip_110_ratios(i); %Scale sides
L_loop(2) =L_loop(2)*Tip_110_ratios(i); %Scale sides

try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(1) = L_trial(1)*Tip_110_ratios(maxindx);
L_trial(2) = L_trial(2)*Tip_110_ratios(maxindx);

L_out = L_trial;

end




%Fits (110) and (001) planes, keeps (111)/(001) ratios the same
function [L_out,translation_out] = dimension_fit_001_pos(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);


[X,Y,L001] = meshgrid(-1:1,-1:1,-1:1); %[X,Y, 001]
step_size = 16;
fit_test = 0;
count_fits =  0;
while count_fits < 10; %10
while fit_test == 0
L_range = [X(:)*pixel_scale_A,Y(:)*pixel_scale_A,L001(:)]*step_size;
for i = 1:27 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
translation_loop = translation_trial;
%L_loop(1) = L_loop(1)*(100+L_range(i,1))/100;
%L_loop(2) = L_loop(2)*(100+L_range(i,1))/100;
L_loop(7) = L_loop(7)*(100+L_range(i,3))/100;
translation_loop = translation_loop + [L_range(i,1),L_range(i,2)];
L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_loop);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
translation_trial = translation_trial+[L_range(maxindx,1),L_range(maxindx,2)];
L_trial(7) = L_trial(7)*(100+L_range(maxindx,3))/100;
L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;
translation_out = translation_trial;

end

function [L_out,translation_out] = dimension_fit_001_pos_v2(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);


[X,Y,L001] = meshgrid(-1:1,-1:1,-1:1); %[X,Y, 001]
step_size = 1;
fit_test = 0;
count_fits =  0;
while count_fits < 10; %10
while fit_test == 0
L_range = [X(:)*pixel_scale_A,Y(:)*pixel_scale_A,L001(:)*2]*step_size;
for i = 1:27 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
translation_loop = translation_trial;
%L_loop(1) = L_loop(1)*(100+L_range(i,1))/100;
%L_loop(2) = L_loop(2)*(100+L_range(i,1))/100;
L_loop(7) = L_loop(7)*(100+L_range(i,3))/100;
translation_loop = translation_loop + [L_range(i,1),L_range(i,2)];
%L_loop = [L_loop(1),L_loop(2),L_loop(3),L_loop(4),L_loop(5),L_loop(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_trial,translation_loop);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
translation_trial = translation_trial+[L_range(maxindx,1),L_range(maxindx,2)];
L_trial(7) = L_trial(7)*(100+L_range(maxindx,3))/100;
%L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.01
    step_size = .01; %Sets lower limit
end
end
L_out = L_trial;
translation_out = translation_trial;

end




%Fits (111) tip ratios
function [ L_out,angles_out] = rotation_fit_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
angles_test = linspace(45,90,181); %angles to test
for j = 1
for i = 1:size(angles_test,2)
L_loop =  L_trial;
angles_loop = angles_trial;
L_projected = max(abs(sqrt(2)*L_trial(1)*cos(angles_test(i)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((angles_test(i)+90)*pi/180))); %Determine projected L
L_loop(1) = L_loop(1)./L_projected*L_loop(1);
L_loop(2) = L_loop(2)./L_projected*L_loop(2);
angles_loop(1) = angles_test(i); %replace with global search
try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_loop,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Ibw(end*3/4:end,:) =0; 
%imshowpair(I,Ibw)
%drawnow
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
    catch
corr_trial(i) = 0;
    end

end
[~,maxindx] = max(corr_trial); %find best guess
L_projected = max(abs(sqrt(2)*L_trial(1)*cos(angles_test(maxindx)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((angles_test(maxindx)+90)*pi/180))); %Determine projected L
L_trial(1) = L_trial(1)./L_projected*L_trial(1);
L_trial(2) = L_trial(1)
angles_trial(1) = angles_test(maxindx);
end
L_out = L_trial;
angles_out = angles_trial;
end


%% 
%Fits (001) tip and rotation
function [L_out,angles_out] = dimension_fit_001_tip_global_rotation(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
Icomp(end/2:end,:) = 0; %mask bottom

linspace(-1,1,41)
%Icomp = imadjust(Icomp);
Tip_001_ratios = linspace(.98,1.02,11); %Ratios to use
angles_test = linspace(45,90,46); %angles to test
[angle_map,tip_map] = meshgrid(angles_test,Tip_001_ratios);
L_ratio = L./L(7); %Compute all ratios compared to (001)
trial_map = [angle_map(:),tip_map(:)];
for i = 1:size(trial_map,1);
angles_loop = angles_trial;
angles_loop(1) = trial_map(i,1); %replace with global search
L_loop =  L_trial;
L_loop(7) =L_trial(7)*trial_map(i,2); %Scale tip
L_projected = max(abs(sqrt(2)*L_trial(1)*cos(trial_map(i,1)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((trial_map(i,1)+90)*pi/180))); %Determine projected L
L_loop(1) = L_loop(1)./L_projected*L_loop(1);
L_loop(2) = L_loop(2)./L_projected*L_loop(2);
L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips

try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_loop,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
Ibw(end/2:end,:) = 0;

corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(7) = L_trial(7)*trial_map(maxindx,2);
angles_trial(1) = trial_map(maxindx,1);

L_projected = max(abs(sqrt(2)*L_trial(1)*cos(trial_map(maxindx,1)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((trial_map(maxindx,1)+90)*pi/180))); %Determine projected L
L_trial(1) = L_trial(1)./L_projected*L_trial(1);
L_trial(2) = L_trial(1)
L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

L_out = L_trial;
angles_out = angles_trial;
end


%% 
%% 
%Fits (001) tip and rotation
function [L_out,angles_out] = dimension_fit_001_tip_global_rotation_fine(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
Icomp(end/2:end,:) = 0;
%Icomp = imadjust(Icomp);
Tip_001_ratios = linspace(.99,1.01,11); %Ratios to use
angles_test = angles_trial(1) + linspace(-5,5,41); %angles to test
[angle_map,tip_map] = meshgrid(angles_test,Tip_001_ratios);
L_ratio = L./L(7); %Compute all ratios compared to (001)
trial_map = [angle_map(:),tip_map(:)];
for i = 1:size(trial_map,1);
angles_loop = angles_trial;
angles_loop(1) = trial_map(i,1); %replace with global search
L_loop =  L_trial;
L_loop(7) =L_trial(7)*trial_map(i,2); %Scale tip
L_projected = max(abs(sqrt(2)*L_trial(1)*cos(trial_map(i,1)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((trial_map(i,1)+90)*pi/180))); %Determine projected L
L_loop(1) = L_loop(1)./L_projected*L_loop(1);
L_loop(2) = L_loop(2)./L_projected*L_loop(2);
L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips

try
[~,shp2d] = nanocrystal_RuO2_shape(L_loop,angles_loop,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
Ibw(end/2:end,:) = 0;

corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(7) = L_trial(7)*trial_map(maxindx,2);
angles_trial(1) = trial_map(maxindx,1);

L_projected = max(abs(sqrt(2)*L_trial(1)*cos(trial_map(maxindx,1)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((trial_map(maxindx,1)+90)*pi/180))); %Determine projected L
L_trial(1) = L_trial(1)./L_projected*L_trial(1);
L_trial(2) = L_trial(1)
L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

L_out = L_trial;
angles_out = angles_trial;
end




