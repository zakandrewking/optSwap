function list = dhRxnList(count)

    if count==31
        list = {'3OAR100';
                '3OAR40';
                '3OAR60';
                '3OAR80';
                'AGPR';
                'AKGDH';
                'ASAD';
                'DHDPRy';
                'EAR40x';
                'EAR60x';
                'EAR80x';
                'G6PDH2r';
                'GAPD';
                'GLUDy';
                'GND';
                'HPYRRx';
                'HSDy';
                'ICDHyr';
                'IPMD';
                'KARA1';
                'KARA2';
                'MDH';
                'MTHFD';
                'NADH16pp';
                'PDH';
                'PGCD';
                'SHK3Dr';
                'TRSARr';
                'ME1';
                'LDH_D';
                'LCARR';};
    elseif count==30
        list = {'3OAR100';
                '3OAR40';
                '3OAR60';
                '3OAR80';
                'AGPR';
                'AKGDH';
                'ASAD';
                'DHDPRy';
                'EAR40x';
                'EAR60x';
                'EAR80x';
                'G6PDH2r';
                'GAPD';
                'GLUDy';
                'GND';
                'HPYRRx';
                'HSDy';
                'ICDHyr';
                'IPMD';
                'KARA1';
                'KARA2';
                'MDH';
                'MTHFD';
                'NADH16pp';
                'PDH';
                'PGCD';
                'SHK3Dr';
        % 'TRSARr'; removed because of loop result
                'ME1';
                'LDH_D';
                'LCARR';};
    elseif count==21
        list = {'ACALD';
                'ALCD2x';
                'ASAD';
                'DHDPRy';
                'FADRx';
                'G6PDH2r';
                'GAPD';
                'GLUDy';
                'GND';
                'HSDy';
                'ICDHyr';
                'IPMD';
                'KARA1';
        % 'MDH'; % removed becuase of loop result, for coupled
                'MTHFD';
                'NADH16pp';
                'PDH';
                'PGCD';
                'SHK3Dr';
                'ME1'; %*
                'LDH_D'; %*
                'LCARR';}; %*
    elseif count==20
        list = {'ACALD';
                'ALCD2x';
                'ASAD';
                'DHDPRy';
                'FADRx';
                'G6PDH2r';
        % 'GAPD';
        % 'GLUDy';% removed becuase of loop result, for yield
                'GND';
                'HSDy';
                'ICDHyr';
                'IPMD';
                'KARA1';
                'MDH'; 
                'MTHFD';
                'NADH16pp';
                'PDH';
                'PGCD';
                'SHK3Dr';
                'ME1'; %*
                'LDH_D'; %*
                'LCARR';}; %*
    else
        error('bad count')
    end

end