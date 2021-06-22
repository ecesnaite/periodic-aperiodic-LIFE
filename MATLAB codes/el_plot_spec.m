function [varargout F] = el_plot_spec_freq_adjust(X,Fs,F_max,varargin)
freqmark = [];

if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for j = 1:2:(length(varargin)-1)
        if ~ischar (varargin{j})
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        switch lower (varargin{j})
            case 'freqmark'
                freqmark = varargin{j+1};  
        end
    end
end
%%
[X_spec,F] = pwelch(X,1*Fs,0,4*Fs,Fs);
X_spec = X_spec(F<=F_max,:);
F = F(F<=F_max);
varargout = X_spec;

hold on
for n = 1:size(X_spec,2)
    command = [ 'disp(''x ' num2str(n) ''')' ];  
    pl(n) = plot(F,10*log10(X_spec(:,n)),'linewidth',1.5, 'ButtonDownFcn', command);  
end
yLimit_spec = get(gca,'YLim');
for k = 1:numel(freqmark)
    hold on, 
    plot(freqmark(k)*[1,1],yLimit_spec,'--','color',rand(1,3));
end

end