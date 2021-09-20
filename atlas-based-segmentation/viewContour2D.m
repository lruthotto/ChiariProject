%% FAIR function for viewing the contour of an image
function varargout = viewContour2D(T,omega,m,color,varargin)
    h  = (omega(2:2:end)-omega(1:2:end))./m;
    xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)';
    [ih, c] = contour(xi(1),xi(2),reshape(T,m)',1, color); axis xy image

    c.LineWidth = 3;

    % the following lines add some nice stuff to the code.
    % if varargin = {'title','FAIR','xlabel','x'}
    % the code evaluates "title('FAIR');xlabel('x');"
    for k=1:2:length(varargin), 
      if ~isempty(varargin{k}), feval(varargin{k},varargin{k+1}); end;
    end;
    if nargout == 1, varargout = {ih};  end;
end