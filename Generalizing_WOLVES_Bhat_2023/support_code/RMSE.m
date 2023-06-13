function rmse = RMSE(y,yhat)
rmse = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error
end
