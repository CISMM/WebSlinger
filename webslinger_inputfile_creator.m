% Webslinger_inputfile_creator
% by Nathan Hudson  07/23
%Updated and modified 07/28/08: NEH

% Open an image of a fibrin network, then you have 3 things to click.
% First, left click on each of the fibrin wall connections.
% Second do the same for the fiber network junctions.  Finally, left click (in
% pairs) on each fiber connection.  For example if wall spot1 and node 3
% are connected by a fiber click on each of those.  To end each clicking
% session, use the right mose button to click on your final point.

function webslinger_inputfile_creator
%NEH: New code for importing image and picking the map to use
[FileName1,PathName1] = uigetfile('*.*','Select an image for Webslinger connections map');
fib_image = imread(strcat(PathName1,FileName1));
fig_location1=[80,80,1100,500]; 
figure('Position',fig_location1);
imshow(fib_image) %NEH: show's image of network
brighten(.6)      %NEH: enhances brightness of image to see fibers.
hold off
%*******************************************************

%*************************************NEH: Webslinger file creation 
%disp('First, pick out the center point') 
%    [x_org,y_org] = ginput(1);

    hold on
disp('Next, click on the spot you pulled from') 
    [x_pull,y_pull] = ginput(1);
    plot(x_pull,y_pull,'c+')
    x_net_min=x_pull(1,1);
    x_net_max=x_pull(1,1);
    y_net_min=x_pull(1,1);
    y_net_max=x_pull(1,1);   
    
%    x_pullN(1,1)=-(x_org(1,1)-x_pull(1,1));
%    y_pullN(1,1)=(y_org(1,1)-y_pull(1,1));
    
hold off

disp('Next, pick out wall connections on figure') 
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
n=0;
hold on
while but == 1
    [xi,yi,but] = ginput(1);
    n=n+1;
    xy_wall(:,n)=[xi,yi];
    plot(xi,yi,'b+')
    text(xi+5,yi+5,num2str(n),'Color',[1,1,0], 'FontSize',8)
%Series of if statements to determine the center of the network    
    if xy_wall(1,n)>=x_net_max 
        x_net_max=xy_wall(1,n);
    end
   if xy_wall(1,n)<=x_net_min 
      x_net_min=xy_wall(1,n);
   end
    if xy_wall(2,n)>=y_net_max
        y_net_max=xy_wall(2,n);
    end
   if xy_wall(2,n)<=y_net_min
      y_net_min=xy_wall(2,n);
   end
    if but ~=1
        large_mass=n;
    end
end
x_org=(x_net_max+x_net_min)/2;
y_org=(y_net_max+y_net_min)/2;
hold off

disp('Next, pick out all network nodes on the figure') 
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
m=0;
hold on
while but == 1
    [xi,yi,but] = ginput(1);
    m=m+1;
    xy_node(:,m)=[xi,yi];
    plot(xi,yi,'g+')
    text(xi+5,yi+5,num2str(m),'Color',[1,0,1], 'FontSize',8)
    small_mass=m;
end
hold off
button=[1;1];
k=1;
disp('Finally, left click all the connected node pairs 2 at a time (INCLUDING THE PULLING POINT)') 
disp('Click the right mouse button on the very last point of the last pair to finish') 
while button(2,1) ~= 3
    [x_pair,y_pair, button] = ginput(2);
%Compare these points to every other point to confirm  
      dist_cut1=100;
      for l=1:large_mass
          dist1=((x_pair(1,1)-xy_wall(1,l))^2+(y_pair(1,1)-xy_wall(2,l))^2)^.5;
           if dist1<= dist_cut1
               junct_mat(k,1)=0;
               junct_mat(k,2)=l;
               dist_cut1=dist1;
           end
      end
      for l=1:small_mass
          dist1=((x_pair(1,1)-xy_node(1,l))^2+(y_pair(1,1)-xy_node(2,l))^2)^.5;
           if dist1<= dist_cut1
               junct_mat(k,1)=1;
               junct_mat(k,2)=l;
               dist_cut1=dist1;
           end
      end 
%PULLED POINT
          dist1=((x_pair(1,1)-x_pull(1,1))^2+(y_pair(1,1)-y_pull(1,1))^2)^.5;
           if dist1<= dist_cut1
               junct_mat(k,1)=1;
               junct_mat(k,2)=0;
               dist_cut1=dist1;
           end
      dist_cut2=100;
      for l=1:large_mass
          dist2=((x_pair(2,1)-xy_wall(1,l))^2+(y_pair(2,1)-xy_wall(2,l))^2)^.5;
           if dist2<= dist_cut2
               junct_mat(k,3)=0;
               junct_mat(k,4)=l;
               dist_cut2=dist2;
           end
      end
      for l=1:small_mass
          dist2=((x_pair(2,1)-xy_node(1,l))^2+(y_pair(2,1)-xy_node(2,l))^2)^.5;
           if dist2<= dist_cut2
               junct_mat(k,3)=1;
               junct_mat(k,4)=l;
               dist_cut2=dist2;
           end
      end
      %PULLED POINT
          dist2=((x_pair(2,1)-x_pull(1,1))^2+(y_pair(2,1)-y_pull(1,1))^2)^.5;
           if dist2<= dist_cut2
               junct_mat(k,1)=1;
               junct_mat(k,2)=0;
               dist_cut2=dist2;
           end
    if button(2,1)==3
        tot_num=k;
    end
    k=k+1;    
end

%Re-locate all points relative to the center of the network and invert x
for i=1:large_mass
    xyW_final(1,i)=-(x_org(1,1)-xy_wall(1,i));
    xyW_final(2,i)=(y_org(1,1)-xy_wall(2,i));
end

for i=1:small_mass
    xyN_final(1,i)=-(x_org(1,1)-xy_node(1,i));
    xyN_final(2,i)=(y_org(1,1)-xy_node(2,i));
end

 x_pullN(1,1)=-(x_org(1,1)-x_pull(1,1)); 
 y_pullN(1,1)=(y_org(1,1)-y_pull(1,1));

outfile = 'loadweb.cfg';

fid = fopen(outfile, 'w+');
fprintf(fid, 'structure { \n');
fprintf(fid, '  mass_damping 3.0 \n \n');
fprintf(fid, '  mass_radius 5.0 \n');

for k=1:large_mass
fprintf(fid, '  mass W%02d  1e100 %03d %03d  0 \n', k, xyW_final(1,k), xyW_final(2,k));
end
fprintf(fid, '\n');
fprintf(fid, '  mass_radius 2.0 \n');
fprintf(fid, '  mass N%02d  1.0 %03d %03d  0 \n', 0, x_pullN(1,1), y_pullN(1,1));
for k=1:small_mass
fprintf(fid, '  mass N%02d  1.0 %03d %03d  0 \n', k, xyN_final(1,k), xyN_final(2,k));
end
fprintf(fid, '\n');
fprintf(fid, '  spring_constant_over_length 10.0 \n');
for k=1:tot_num
    node_check1=junct_mat(k,1);
    node_check2=junct_mat(k,3);
    if (node_check1==0) && (node_check2==0)
             fprintf(fid, '  spring W%02d W%02d \n', junct_mat(k,2), junct_mat(k,4));
    elseif (node_check1==0) && (node_check2==1)
             fprintf(fid, '  spring W%02d N%02d \n', junct_mat(k,2), junct_mat(k,4));
    elseif (node_check1==1) && (node_check2==0)
             fprintf(fid, '  spring N%02d W%02d \n', junct_mat(k,2), junct_mat(k,4));
    else
             fprintf(fid, '  spring N%02d N%02d \n', junct_mat(k,2), junct_mat(k,4));
    end
end
fprintf(fid, '} \n');
stat=fclose(fid);
%*************************************NEH: Webslinger file creation 

end