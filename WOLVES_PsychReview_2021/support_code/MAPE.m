function mape = MAPE (y, yhat)
    mape  = mean(abs((y - yhat)./y))*100; % Mean Absolute Percentage Error
end