
function[] = multilayer_film_plot_RTA(lambda, theta, R, T, A, pol)

if (length(theta)==1)
    for p = 1:length(pol)
        figure;
        semilogx(lambda, [R{p} T{p} A{p}]);
        xlim([min(lambda) max(lambda)]);
        titleStr = sprintf('%s polarization, angle = %g degrees', pol{p}, theta);
        title(titleStr);
        xlabel('Wavelength');
        ylabel('Refl / Trans / Abs');
        legend({'R','T','A'});
    end
else
    for p = 1:length(pol)
        figure;
        imagesc(theta*180/pi, lambda, R{p});
        titleStr = sprintf('Reflection, %s Polarization', pol{p});
        title(titleStr);
        xlabel('Angle (degrees)')
        ylabel('Wavelength');
        
        figure;
        imagesc(theta*180/pi, lambda, T{p});
        titleStr = sprintf('Transmission, %s Polarization', pol{p});
        title(titleStr);
        xlabel('Angle (degrees)')
        ylabel('Wavelength');
        
        figure;
        imagesc(theta*180/pi, lambda, A{p});
        titleStr = sprintf('Absorption, %s Polarization', pol{p});
        title(titleStr);
        xlabel('Angle (degrees)')
        ylabel('Wavelength');
    end
end

end