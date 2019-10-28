clc;
close all;
clear all;

I_frame = imread("imageA.bmp");
P_frame = imread("imageB.bmp");
tic%��ʱ
K =1;%ѹ������

%��׼������
Z = [16,11,10,16,24,40,51,61
    12,12,14,19,26,58,60,55
    14,13,16,24,40,57,69,56
    14,17,22,29,51,87,80,62
    18,22,37,56,68,109,103,77
    24,35,55,64,81,104,113,92
    49,64,78,87,103,121,120,101
    72,92,95,98,112,100,103,99];
R=8; %������Χ ������ǶԿ����������ĸ����� ��8�����ص�ʮ����������п�ƥ��


%**************************************
%�����ȶ�i֡����ѹ����

i = I_frame;%��ȡͼ��
I = i;
[h0,w0]=size(i);
i =padding(i);%padding 

%����������ͼ��ѹ������
[H1,W1]=size(i);
nH=H1/8;%������ٿ�
nW=W1/8;%������ٿ�
idx=0;
%block(1:8,1:8,1:nH*nW)=zeros(8,8,nH*nW);%Ԥ�ȷ���ռ䣬���Լӿ������ٶ�
%Qblock(1:8,1:8,1:nH*nW)=zeros(8,8,nH*nW);
Zline1 = [];%Ԥ����ռ� �洢������
for j=1:nH
    for k=1:nW
    idx=idx+1;
    
    block(:,:,idx)=i(((j-1)*8+1:(j-1)*8+8),((k-1)*8+1:(k-1)*8+8));%�ֿ�
    DCTblock(:,:,idx) = dct2(block(:,:,idx));%����dct�任
    Qblock(:,:,idx) =round(DCTblock(:,:,idx)./(K*Z));%����������֮����./ 
    
    n=0;%���һ������λ��
    m=64;%����
    Zline=zigzag(Qblock(:,:,idx));%zigzag %ȥ������ȫ�����ֵ
        while m~=0
            if Zline(m)~=0
                n=m;
                break;
            end 
        m=m-1;
        end
    Zline1 =[Zline1,n,Zline(1:n)];
    end
end      
%���ˣ���Ҫѹ�������Ѿ���� 
%�ļ�ͷ
%header 40λ
%����
%Height 16λ
%Width  16λ
%ѹ����8λ
header =[h0,w0,K,nH,nW,0,0,0];
Seq =[header,double(Zline1)];
%ѹ�����

Seq0 = Seq; %������������
%��ȡheader 
rH = Seq0(4);
rW = Seq0(5);
rh0 = Seq0(1); %ԭʼ��Сh
rw0 = Seq0(2); %ԭʼ��Сw
rK = Seq0(3);  %ѹ����
idx=0;
%�ļ�ͷ��ȡ���֮�� ���Ǵ�����ܵı�����
base = 9;
for j = 1:rH
    for k =1:rW
        idx = idx +1;
        ra =Seq0(base+1:base+Seq0(base));
        rQblock(:,:,idx)=izigzag(ra);
        rDCTblock(:,:,idx)=rQblock(:,:,idx).*(rK*Z); %�����Դ֮��
        rblock(:,:,idx)=idct2(rDCTblock(:,:,idx));
        rimage((j-1)*8+1:j*8,(k-1)*8+1:k*8)=rblock(:,:,idx);
        base =Seq0(base)+base+1;
    end
end
rimage =round(rimage);
rimage_out = rimage(1:rh0,1:rw0);
rimage_out = uint8(rimage_out);
figure;
subplot(121),imshow(I),title("the image"); %�۲�һ��ԭͼ
subplot(122),imshow(rimage_out),title("image after being compressed");

%********************************************

%һЩ��������
origin_image = h0*w0*8;%bit
compress_image = length(Seq)*8;%bit ���Թ���
CR = origin_image./compress_image;
PSNR=10*log10(255*255/mean(mean((double(I)-double(rimage)).^2)));

disp(['When K=',num2str(K),':']);
disp(['for I_frame : ']);
disp(['PSNR                 ',num2str(PSNR)]);
disp(['Original Bit:         ',num2str(origin_image),' bit']);
disp(['Compressed Bit:       ',num2str(compress_image),' bit']);
disp(['Compression Ratios for I_frame:   ',num2str(CR)]);
%**************************************


% p֡��������
iP = padding(P_frame); 
[Hp,Wp] = size(iP);
nhp = Hp/8;
nwp = Wp/8;

ZlineP =[]; %used to keep the compressed seq
cellvector=cell(nhp,nwp); % used to keep the vector  for each block
i =padarray(i,[R,R],'replicate','both'); % I֡��Χ������P֡�߽�����
rimage = padarray(rimage,[R,R],'replicate','both'); 
for x = 1:nhp
    for y = 1:nwp
        pblock = iP((1+8*(x-1):8*x),(1+8*(y-1):8*y));
        %ƥ��
        vector = Bestmatch(pblock, i, x, y, R ); %ssD
        %vector = Bestmatch(pblock, I_frame, x, y, R ,SAD);
        cellvector{x,y} = vector;
        r = vector(1,1);
        c = vector(1,2);
        mblock =i((1+8*(x-1)+R+r:8*x+R+r),(1+8*(y-1)+R+c:8*y+R+c)); % �ҵ������ƥ���
        Resblock = pblock-mblock;
        DCTResblock=dct2(Resblock);
        QDCTResblock =round( DCTResblock./(K*Z));
        %zigzag
        ZlineP1 = zigzag(QDCTResblock);
        
        n=0;%���һ������λ��
        m=64;%����
        while m~=0
            if ZlineP1(m)~=0
                n=m;
                break;
            end 
            m=m-1;
        end
        ZlineP =[ZlineP,n,ZlineP1(1:n)];%ѹ���������
        
        % ��������˵����ʵ��Ӧ���У�������Ҫ����ļ�ͷ�Ա��ѹ��just like what we do in jpeg 
        % Ϊ�˼��ٹ������������ѭ���ṹ�У����Ǿͽ��н�ѹ��
        IQDCTResblock = round(QDCTResblock.*(K*Z));
        IDCTResblock = idct2(IQDCTResblock);
        recoverblock = double(IDCTResblock) + double(rimage((1+8*(x-1)+R+r:8*x+R+r),(1+8*(y-1)+R+c:8*y+R+c)));
        rP_frame((1+8*(x-1):8*x),(1+8*(y-1):8*y)) = recoverblock;
    end
end
toc
%����ѹ����
ImageB=Hp*Wp*8;
Pframe=length(ZlineP)*8+nhp*nwp*2*2*log2(R); %���е�ÿһ��Ԫ�ض�Ϊ8bit��Ԫ����ÿ��Block���˶�ʸ��Ϊ2*log2(R)bit
CR1=ImageB/Pframe;
PSNR1=10*log10(255*255/mean(mean((double(iP)-double(rP_frame)).^2))); %P֡��ֵ�����
disp(['When K=',num2str(K),':']);
disp(['for P_frame : ']);
disp(['PSNR:                 ',num2str(PSNR1)]);
disp(['ImageB Bit:           ',num2str(ImageB),' bit']);
disp(['PFrame Bit:           ',num2str(Pframe),' bit']);
disp(['Compression Ratios:   ',num2str(CR1)]);

figure,
imshow(uint8(rP_frame)),title("ImageB Ԥ��ѹ����ѹ��");


function y=zigzag(a)
    zz=[1,2,9,17,10,3,4,11,18,25,33,26,19,12,5,6,13,20,27,34,41,49,42,35,28,21,14,7,...
        8,15,22,29,36,43,50,57,58,51,44,37,30,23,16,24,31,38,45,52,59,60,53,46,39,32,...
        40,47,54,61,62,55,48,56,63,64];
    aa = reshape(a,1,64);%ת����������
    y=aa(zz);% �������еĵ�zz���Ƕ�Ӧ��y�е�
end

function y =izigzag(a)
        zz=[1,2,9,17,10,3,4,11,18,25,33,26,19,12,5,6,13,20,27,34,41,49,42,35,28,21,14,7,...
        8,15,22,29,36,43,50,57,58,51,44,37,30,23,16,24,31,38,45,52,59,60,53,46,39,32,...
        40,47,54,61,62,55,48,56,63,64];
        b = zeros(1,64); %���ٴ洢�ռ�
        for i=1:length(a)
            b(zz(i))=a(i);
        end
        y =reshape(b,8,8);
end

function error = SSD(blockA , blockB)
    cost = double(double(blockA) -double(blockB)).*double(double(blockA)-double(blockB));
    error =sum(cost(:));
end

function error = SAD(blockA , blockB)
    cost =abs(blockA -blockB);
    error =sum(cost);
end

function image = padding(i)

    [h0,w0]=size(i);%��ȡ�ߴ���Ϣ
    %padding
    H = ceil(h0/8)*8;
    W = ceil(w0/8)*8;
    %padding h using the last column
    if mod(h0,8)~=0
        for m=h0:H
            i(m,:)=i(h0,:);
        end
    end
 %padding w using the last line
 
    if mod(w0,8)~=0
        for m=w0:W
            i(:,m)=i(:,w0);
        end
    end
    image =i;
%pading ���
end

function [vector] =Bestmatch(pblock, I, row, col, R)
    %R��������Χ
    vector =[0,0]; %�Ե�ǰ��Ϊ����ԭ�㡣
    %��Ѱ�ҵ���ͬλ�õĿ顣
   % I=padarray(I,[R,R],'replicate','both'); % I֡��Χ������P֡�߽�����
    Iblock=I((1+8*(row-1):8*row),(1+8*(col-1):8*col));
    % Ŀ�� �ҵ���С��SAD or SSD
    %SSD�汾
        error = SSD(pblock,Iblock); 
        %error = SAD(pblock,Iblock);   
    for i=-R+1:R-1
        for j= -R+1:R-1
                Iblock=I((1+8*(row-1)+i+R:8*row+R+i),(1+8*(col-1)+R+j:8*col+R+j));
                errorNow = SSD(pblock,Iblock);
                if errorNow < error
                    vector = [i,j];
                    error = errorNow;
                end 
        end
    end
           
end