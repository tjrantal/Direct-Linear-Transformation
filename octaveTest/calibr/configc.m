function sys=configc(name)
%CONFIGC Loads system configuration information based on given name. 
%User may edit here new setup data.
if strcmp(name,'sony')
  sys = [
      768,     %number of pixels in horizontal direction
      576,     %number of pixels in vertical direction
      6.2031,  %effective CCD chip size in horizontal direction
      4.6515,  %effective CCD chip size in vertical direction
      8.5,     %nominal focal length
      5,       %radius of the circular control points
      0,       %for future expansions
      0,
      0,
      abs(name)'
  ];
  return;
end

if strcmp(name,'sonyz')
  sys = [
      768,     %number of pixels in horizontal direction
      576,     %number of pixels in vertical direction
      6.2031,  %effective CCD chip size in horizontal direction
      4.6515,  %effective CCD chip size in vertical direction
      8.5,     %nominal focal length
      0,       %radius of the circular control points
      0,       %for future expansions
      0,
      0,
      abs(name)'
  ];
  return;
end
      
if strcmp(name,'pulnix')
  sys = [
      512,     %number of pixels in horizontal direction
      512,     %number of pixels in vertical direction
      4.3569,  %effective CCD chip size in horizontal direction
      4.2496,  %effective CCD chip size in vertical direction
      16,      %nominal focal length
      15,      %radius of the circular control points
      0,       %for future expansions
      0,
      0,
      abs(name)'
  ];
  return;
end

error('Unknown camera type')
