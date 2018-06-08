close all; clear; clc;

node = load('Coord.dat');
elem = load('Connectivity.dat');
%elem = elem + (1 + 0*elem);
num_nodes = size(node,1);
num_elem = size(elem,1);
U = load('U.dat');
U = U';

uu(num_nodes,1)=0;
vv(num_nodes,1)=0;
for i=1:num_nodes
    uu(i,1)=U(2*i-1,1);
    vv(i,1)=U(2*i,1);
end

node_changed=node+[uu,vv];
figure
for i=1:num_elem
    for j=1:4
        pp(j,:)=node_changed(elem(i,j),:);
    end
    patch(pp(:,1),pp(:,2),[1 0 0])
    hold on    
    grid on
end