function plotImageGray(img,str0)
        figure;
        imagesc(db(abs(img)/max(abs(img(:)))),[-60 0]); % colormap('gray')
        title(str0)
        xlabel('Azimuth Samples');ylabel('Range Samples');
        axis xy;
end
