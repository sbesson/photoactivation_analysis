function originalFile = updateOriginalFile(session, originalFile, filePath, varargin)
% UPDATEORIGINALFILE Replace a file with new content on the OMERO server
%
%    originalFile = updateOriginalFile(session, originalFile, filePath)
%    replaces the content of the original file specified by filePath
%
%    originalFile = updateOriginalFile(session, originalFile, filePath,
%    nTries) also specifies the number of tries of the upload fails.
%
%    Examples:
%
%        originalFile = updateOriginalFile(session, originalFile, filePath)
%        originalFile = updateOriginalFile(session, originalFile, filePath, nTries)
%
% See also: createFileAnnotation

% Copyright (C) 2013 University of Dundee & Open Microscopy Environment.
% All rights reserved.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

% Input check
ip = inputParser;
ip.addRequired('session');
ip.addRequired('originalFile', @(x) isa(x, 'omero.model.OriginalFile'));
ip.addRequired('filePath', @(x) exist(x, 'file') == 2);
ip.addOptional('nTries', 1, @isscalar);
ip.parse(session, originalFile, filePath, varargin{:});

% Read absolute input file path
[~, f] = fileattrib(filePath);
absolutePath = f.Name;

% Compare checksums of file client-side verus server-side
[originalFile, hasher] = updateFile(session, originalFile, absolutePath);
clientHash = char(hasher.checksumAsString());
serverHash = char(originalFile.getSha1().getValue());
iTry = 1;

% Allow to re-upload of data if mismatching checksums
msg = 'File checksum mismatch on upload: %s (client has %s, server has %s)';
while iTry < ip.Results.nTries && ~strcmp(clientHash, serverHash)
    fprintf(1, [msg '. Retrying...\n'], filePath, clientHash, serverHash);
    iTry = iTry + 1;
    [originalFile, hasher] = updateFile(session, originalFile, absolutePath);
    clientHash = char(hasher.checksumAsString());
    serverHash = char(originalFile.getSha1().getValue());
end

% Check the file has been correctly uploaded
assert(isequal(clientHash, serverHash), msg, filePath, clientHash, serverHash);

function [originalFile, hasher] = updateFile(session, originalFile, absolutePath)

% Create java io File
[path, name, ext] = fileparts(absolutePath);
fileLength = length(java.io.File(absolutePath));

% Create original file
originalFile.setName(rstring([name ext]));
originalFile.setPath(rstring(path));
originalFile.setSize(rlong(fileLength));

% now we save the originalFile object
updateService = session.getUpdateService();
originalFile = updateService.saveAndReturnObject(originalFile);

% Initialize provider to compute client-side checksum
checksumProviderFactory = ome.util.checksum.ChecksumProviderFactoryImpl;
sha1= ome.util.checksum.ChecksumType.SHA1;
hasher = checksumProviderFactory.getProvider(sha1);

% Initialize the service to load the raw data
rawFileStore = session.createRawFileStore();
rawFileStore.setFileId(originalFile.getId().getValue());

%code for small file.
fid = fopen(absolutePath);
byteArray = fread(fid,[1, fileLength], 'uint8');
rawFileStore.write(byteArray, 0, fileLength);
hasher.putBytes(byteArray);
fclose(fid);

% Save and close the service
originalFile = rawFileStore.save();
rawFileStore.close();