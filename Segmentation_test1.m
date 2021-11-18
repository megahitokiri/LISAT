%%
cd Tray17A_SOS20E02-002_20211029_14341386
map = jet(numberOfImagesToProcess); 
%colormap(map);
%%
%%
%L* SPACE START...
figure();
for k = 1 : numberOfImagesToProcess
	% Read in this one file.
	baseFileName = imageFiles(k).name;
	fullFileName = fullfile(imagesFolder, baseFileName);
	rgbImage = imread(fullFileName);
   
    % Program Definition
    leaf = rgbImage;
    %Transform into lab space
    lab_leaf = rgb2lab(leaf);
        l = lab_leaf(:,:,1);
        a = lab_leaf(:,:,2);
        b = lab_leaf(:,:,3);

    %Compose histogram results for l
    [counts,edges] = histcounts(l,100);
         edges= edges(1:100);
         l_results = [edges;counts].';  % vertical concatenation & transpose
    
    %Applies a Savitzky-Golay finite impulse response (FIR) smoothing filter of polynomial order
    order = 3; framelen = 13;
    l_sgf = sgolayfilt(l_results,order,framelen);
    l_PeakPos=l_sgf(:,1);
    l_PeakSig=l_sgf(:,2);
    
    %Plotting Smooth transformation
    plot(l_PeakPos,l_PeakSig,"Color",map(k,:),'DisplayName',baseFileName);
    title('L* Space Savitzky-Golay smoothed (Superimposed)', 'FontSize', fontSize);
    alpha(0.5);
    hold on
    
    [l_pks,l_locs] = findpeaks(l_PeakSig,l_PeakPos);
    l_pk_locs = [l_locs,l_pks];
    
    %selecting local minimums and bimodal based on Whithe-Thresholding
end    
hold off
legend

%%
%A* SPACE START...
figure();
for k = 1 : numberOfImagesToProcess
	% Read in this one file.
	baseFileName = imageFiles(k).name;
	fullFileName = fullfile(imagesFolder, baseFileName);
	rgbImage = imread(fullFileName);
   
    % Program Definition
    leaf = rgbImage;
    %Transform into lab space
    lab_leaf = rgb2lab(leaf);
        l = lab_leaf(:,:,1);
        a = lab_leaf(:,:,2);
        b = lab_leaf(:,:,3);

    %Compose histogram results for a
    [counts,edges] = histcounts(a,100);
         edges= edges(1:100);
         a_results = [edges;counts].';  % vertical concatenation & transpose
    
    %Applies a Savitzky-Golay finite impulse response (FIR) smoothing filter of polynomial order
    order = 3; framelen = 13;
    a_sgf = sgolayfilt(a_results,order,framelen);
    a_PeakPos=a_sgf(:,1);
    a_PeakSig=a_sgf(:,2);
    
    %Plotting Smooth transformation
    plot(a_PeakPos,a_PeakSig,"Color",map(k,:),'DisplayName',baseFileName);
    title('A Space Savitzky-Golay smoothed (Superimposed)', 'FontSize', fontSize);
    alpha(0.5);
    hold on
    
    [a_pks,a_locs] = findpeaks(a_PeakSig,a_PeakPos);
    a_pk_locs = [a_locs,a_pks];
    
end    
hold off
legend
%%
%B* SPACE START...
figure();
for k = 1 : numberOfImagesToProcess
	% Read in this one file.
	baseFileName = imageFiles(k).name;
	fullFileName = fullfile(imagesFolder, baseFileName);
	rgbImage = imread(fullFileName);
   
    % Program Definition
    leaf = rgbImage;
    %Transform into lab space
    lab_leaf = rgb2lab(leaf);
        l = lab_leaf(:,:,1);
        a = lab_leaf(:,:,2);
        b = lab_leaf(:,:,3);
    %Compose histogram results for b
    [counts,edges] = histcounts(b,100);
         edges= edges(1:100);
         b_results = [edges;counts].';  % vertical concatenation & transpose
    
    %Applies a Savitzky-Golay finite impulse response (FIR) smoothing filter of polynomial order
    order = 3; framelen = 13;
    b_sgf = sgolayfilt(b_results,order,framelen);
    b_PeakPos=b_sgf(:,1);
    b_PeakSig=b_sgf(:,2);
    
    plot(b_PeakPos,b_PeakSig,"Color",map(k,:),'DisplayName',baseFileName);
    title('B* Space Savitzky-Golay smoothed (Superimposed)', 'FontSize', fontSize);
    alpha(0.5);
    hold on
    
    [b_pks,b_locs] = findpeaks(b_PeakSig,b_PeakPos);
    b_pk_locs = [b_locs,b_pks];
    
end    
hold off
%%
%Individual Analysis
for k = 1 : numberOfImagesToProcess
    figure();
	% Read in this one file.
	baseFileName = imageFiles(k).name;
	fullFileName = fullfile(imagesFolder, baseFileName);
	rgbImage = imread(fullFileName);
   
    % Program Definition
    leaf = rgbImage;
    %Transform into lab space
    lab_leaf = rgb2lab(leaf);
        l = lab_leaf(:,:,1);
        a = lab_leaf(:,:,2);
        b = lab_leaf(:,:,3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compose histogram results for L*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2, 3, 1);
        %Compose histogram results for L
        [counts,edges] = histcounts(l,100);
             edges= edges(1:100);
             l_results = [edges;counts].';  % vertical concatenation & transpose
        
        %Applies a Savitzky-Golay finite impulse response (FIR) smoothing filter of polynomial order
        order = 3; framelen = 13;
        l_sgf = sgolayfilt(l_results,order,framelen);
        l_PeakPos=l_sgf(:,1);
        l_PeakSig=l_sgf(:,2);
        
        %Plotting Smooth transformation
        plot(l_PeakPos,l_PeakSig,"Color",map(k,:),'DisplayName',baseFileName);
        title(strcat('L-Space:',baseFileName), 'FontSize', fontSize);    
        [l_pks,l_locs,l_widths,l_proms] = findpeaks(l_PeakSig,l_PeakPos,'MinPeakProminence',50000,'Annotate','extents','MinPeakWidth',2);
        l_pk_locs = [l_locs,l_pks,l_widths,l_proms];
        
        findpeaks(l_PeakSig,l_PeakPos,'MinPeakProminence',50000,'Annotate','extents','MinPeakWidth',2);
        title(strcat('L-Space:',baseFileName), 'FontSize', fontSize);    
    
        text(l_locs+.02,l_pks,num2str((1:numel(l_pks))'))
    
        l_pk_locs_table = array2table(l_pk_locs);
        %selecting local minimums and bimodal based on Whithe-Thresholding
    	legend

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compose histogram results for A*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2, 3, 2);
        %Compose histogram results for a
        [counts,edges] = histcounts(a,100);
             edges= edges(1:100);
             a_results = [edges;counts].';  % vertical concatenation & transpose
        
        %Applies a Savitzky-Golay finite impulse response (FIR) smoothing filter of polynomial order
        order = 3; framelen = 13;
        a_sgf = sgolayfilt(a_results,order,framelen);
        a_PeakPos=a_sgf(:,1);
        a_PeakSig=a_sgf(:,2);
        
        %Plotting Smooth transformation
        plot(a_PeakPos,a_PeakSig,"Color",map(k,:),'DisplayName',baseFileName);
        title(strcat('A-Space:',baseFileName), 'FontSize', fontSize);    
        [a_pks,a_locs,a_widths,a_proms] = findpeaks(a_PeakSig,a_PeakPos,'MinPeakProminence',50000,'Annotate','extents','MinPeakWidth',2);
        a_pk_locs = [a_locs,a_pks,a_widths,a_proms];
        
        findpeaks(a_PeakSig,a_PeakPos,'MinPeakProminence',5000,'Annotate','extents','MinPeakWidth',2);
        title(strcat('A-Space:',baseFileName), 'FontSize', fontSize);    
    
        text(a_locs+.02,a_pks,num2str((1:numel(a_pks))'))
    
        a_pk_locs_table = array2table(a_pk_locs);
        %selecting local minimums and bimodal based on Whithe-Thresholding
        legend

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Compose histogram results for B*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(2, 3, 3);
        %Compose histogram results for b
        [counts,edges] = histcounts(b,100);
             edges= edges(1:100);
             b_results = [edges;counts].';  % vertical concatenation & transpose
        
        %Applies a Savitzky-Golay finite impulse response (FIR) smoothing filter of polynomial order
        order = 3; framelen = 13;
        b_sgf = sgolayfilt(b_results,order,framelen);
        b_PeakPos=b_sgf(:,1);
        b_PeakSig=b_sgf(:,2);
        
        %Plotting Smooth transformation
        plot(b_PeakPos,b_PeakSig,"Color",map(k,:),'DisplayName',baseFileName);
        title(strcat('B-Space:',baseFileName), 'FontSize', fontSize);    
        [b_pks,b_locs,b_widths,b_proms] = findpeaks(b_PeakSig,b_PeakPos,'MinPeakProminence',50000,'Annotate','extents','MinPeakWidth',2);
        b_pk_locs = [b_locs,b_pks,b_widths,b_proms];
        
        findpeaks(b_PeakSig,b_PeakPos,'MinPeakProminence',50000,'Annotate','extents','MinPeakWidth',2);
        title(strcat('B-Space:',baseFileName), 'FontSize', fontSize);    
    
        text(b_locs+.02,b_pks,num2str((1:numel(b_pks))'))
    
        b_pk_locs_table = array2table(b_pk_locs);
        sorted_b_pk_locs_table = sortrows(b_pk_locs_table, 'b_pk_locs4','descend'); % sort the table by 'prominence Size'
        filter_b_pk_locs_table = sorted_b_pk_locs_table(1:2,:);
        b_pk_locs_table = sortrows(filter_b_pk_locs_table, 'b_pk_locs1','ascend');
        %selecting local minimums and bimodal based on Whithe-Thresholding
        legend

end    

%%
cd ..
imsubstract


