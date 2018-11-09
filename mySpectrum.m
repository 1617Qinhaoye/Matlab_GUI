function  result=mySpectrum(f)
    result=uint8(255*mat2gray(log(abs(fftshift(fft2((f))))+1)));%进行傅里叶变换,并中心化,对其值进行log变换,以便于显示
end