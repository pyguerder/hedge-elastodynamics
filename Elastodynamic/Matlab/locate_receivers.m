function [nodes, elements] = locate_receivers(receiver_x, receiver_y, ...
                                              space_x, space_y)
% Find node numbers (row index in x) and element numbers (column index in
% x) of the nodes close to the receiver.

% After this, we also have to find in the corresponding reference triangle its
% coordinates (r0,s0) (see page 81).  Note that although in most cases this
% function can return more than one node, in Elastics2D we only take the first
% node to extract the signal from.

global x y

[nodes, elements] = find(  x <= receiver_x + space_x ...
                         & x >= receiver_x - space_x ...
                         & y <= receiver_y + space_y ...
                         & y >= receiver_y - space_y );

%Pierre-Yves : essayer nonzero() en Python