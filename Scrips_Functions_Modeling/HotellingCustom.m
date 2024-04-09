function [T_squared, F_stats, df1, df2, p_value_determined] = HotellingCustom(dataToClassifyAllButOthers_Paroxysmal, dataToClassifyAllButOthers_Persistent)

cov_Paroxysmal = cov(dataToClassifyAllButOthers_Paroxysmal);
cov_Persistent = cov(dataToClassifyAllButOthers_Persistent);


S_both = ( ((size(dataToClassifyAllButOthers_Paroxysmal,1)-1).*(cov_Paroxysmal)) + ((size(dataToClassifyAllButOthers_Persistent,1)-1).*(cov_Persistent)) )./( size(dataToClassifyAllButOthers_Paroxysmal,1) + size(dataToClassifyAllButOthers_Persistent,1) - 2 );
invS_both = inv(S_both);

MD_both = sqrt( ( mean(dataToClassifyAllButOthers_Paroxysmal,1) - mean(dataToClassifyAllButOthers_Persistent,1) ) * invS_both * ( mean(dataToClassifyAllButOthers_Paroxysmal,1) - mean(dataToClassifyAllButOthers_Persistent,1) )' );

T_squared = (( size(dataToClassifyAllButOthers_Persistent,1) * size(dataToClassifyAllButOthers_Paroxysmal,1) )./( size(dataToClassifyAllButOthers_Persistent,1) + size(dataToClassifyAllButOthers_Paroxysmal,1) ))*MD_both.^2;

F_stats = ( ( size(dataToClassifyAllButOthers_Persistent,1) + size(dataToClassifyAllButOthers_Paroxysmal,1) - size(dataToClassifyAllButOthers_Persistent,2) - 1 )./( size(dataToClassifyAllButOthers_Persistent,2)*( size(dataToClassifyAllButOthers_Persistent,1) + size(dataToClassifyAllButOthers_Paroxysmal,1) - 2 ) ) ).*( T_squared );

df1 = size(dataToClassifyAllButOthers_Persistent,2);
df2 = size(dataToClassifyAllButOthers_Persistent,1) + size(dataToClassifyAllButOthers_Paroxysmal,1) - size(dataToClassifyAllButOthers_Persistent,2) - 1;

x_stat = 0:0.01:F_stats*2;
y_stat = fpdf(x_stat,df1,df2);
figure;
plot(x_stat,y_stat);

indInterestF = find( abs(x_stat - F_stats) <= .01 );
indInterestF = indInterestF(1);
p_value_determined = trapz(x_stat(indInterestF:end),y_stat(indInterestF:end));

[T_squared, F_stats, df1, df2, p_value_determined]

end


