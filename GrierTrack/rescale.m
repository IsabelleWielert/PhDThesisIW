function [t_out, v_out, v_sm_out] = rescale(t,v,v_sm,delta_t)

t = t - delta_t;
ind = find(t>=0);
t_out = t(ind);
v_out = v(ind);
v_sm_out = v_sm(ind);
