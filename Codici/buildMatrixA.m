function A = buildMatrixA(nComponents)
% Builds matrix A required to solve the linear system Ax = b where:

    
    % --- Creazione della matrice A (inizializzata a zero)
    A = zeros(3 * nComponents, 3 * nComponents);
    
    % --- Completamento della matrice A con le condizioni di equilibrio
    for i = 1:nComponents
        % Equilibrio assiale (forze lungo x)
        A(i, i) = 1; % Forza assiale del componente i
        if i < nComponents
            A(i, i + 1) = -1; % Trasmissione della forza assiale al componente successivo
        end
        
        % Equilibrio trasversale (forze lungo z)
        A(nComponents + i, nComponents + i) = 1; % Forza trasversale del componente i
        if i < nComponents
            A(nComponents + i, nComponents + i + 1) = -1; % Trasmissione della forza trasversale
        end
        
        % Equilibrio momenti (momenti torcenti)
        A(2 * nComponents + i, 2 * nComponents + i) = 1; % Momento torcenti del componente i
        if i < nComponents
            A(2 * nComponents + i, 2 * nComponents + i + 1) = -1; % Trasmissione del momento
        end
    end
    
    % --- Applichiamo le condizioni di continuità per le forze e i momenti
    for i = 1:nComponents - 1
        % Continuity in direzione assiale
        A(i, i + 1) = -1;
        A(i + 1, i) = 1;
        
        % Continuity in direzione trasversale
        A(nComponents + i, nComponents + i + 1) = -1;
        A(nComponents + i + 1, nComponents + i) = 1;
        
        % Continuity nei momenti
        A(2 * nComponents + i, 2 * nComponents + i + 1) = -1;
        A(2 * nComponents + i + 1, 2 * nComponents + i) = 1;
    end
end
