function tw = setupTwitty

    credentials.ConsumerKey = 'M8AxYB091hiaQ7HTHZgtA';
    credentials.ConsumerSecret = 'his0u2LEmXO01RU4MTqXY1zbmJRP4PakfJulYNveA';
    credentials.AccessToken = '735579048-Iijijhb9LbOGi1l55VOK4Z0ylsgxjrelZf39tac0';
    credentials.AccessTokenSecret = 'U0kV7zuN7UJdUcgN4pVGgc4jya7AmS7spNwl6PR2o';
    
    try
        tw = twitty(credentials);
    catch
        tw = false;
    end
end