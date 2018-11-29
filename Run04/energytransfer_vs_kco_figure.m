h1 = openfig('I1Pdf8.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
h2 = openfig('I1Pdf16.fig','reuse'); % open figure
ax2 = gca; % get handle to axes of figure
h3 = openfig('I1Pdf32.fig','reuse'); % open figure
ax3 = gca; % get handle to axes of figure
h4 = openfig('I1Pdf64.fig','reuse'); % open figure
ax4 = gca; % get handle to axes of figure
h5 = openfig('I1Pdf128.fig','reuse'); % open figure
ax5 = gca; % get handle to axes of figure
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots
h6 = figure; %create new figure
s1 = subplot(2,3,1); %create and get handle to the subplot axes
s2 = subplot(2,3,2);
s3 = subplot(2,3,3);
s4 = subplot(2,3,4);
s5 = subplot(2,3,5);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
fig3 = get(ax3,'children');
fig4 = get(ax4,'children');
fig5 = get(ax5,'children');
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
copyobj(fig3,s3);
copyobj(fig4,s4);
copyobj(fig5,s5);