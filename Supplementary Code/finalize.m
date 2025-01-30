function finalize(ofo)
h = ancestor(ofo, 'figure');
% Reset the mouse scroll wheel callback.
h.WindowScrollWheelFcn = [];
% Save finalized set of points.
pos = ofo.Position;
% Delete and create a new Freehand ROI with the new |Position| value.
delete(ofo);
drawfreehand(gca, 'Position', pos);
end