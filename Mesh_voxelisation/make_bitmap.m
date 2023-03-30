function [G_S_bits, G_H_bits] = make_bitmap(G_greyscale, rho_e, G_S, G_H, U_xy, U_z)

G_H_bits = repelem(G_greyscale,U_xy,U_xy,U_z);

for rhoi = 1 : length(rho_e)
    rho = rho_e(rhoi);
    rho(rho==0.001) = 0;
    inds = find( G_H_bits == rho );
    inds1 = randperm(length(inds), round(length(inds)*rho));
    G_H_bits(inds) = 0;
    G_H_bits(inds(inds1)) = 1;
end
G_H_bits = logical(round(G_H_bits));
G_S_bits = ((~G_H_bits));

G_H_bits = G_H_bits&(repelem(G_H+G_S,U_xy,U_xy,U_z));
G_S_bits = G_S_bits&(repelem(G_H+G_S,U_xy,U_xy,U_z));