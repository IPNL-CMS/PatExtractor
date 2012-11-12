{

  int total = tree->Draw("mtt", "");
  int kfConverged = tree->Draw("kf_converged", "kf_converged");
  int associable = tree->Draw("associable", "associable && kf_converged");
  int goodCombinaison = tree->Draw("good_combinaison", "kf_converged && good_combinaison && associable");
  int wellPlaced = tree->Draw("combinaison_well_placed_after_kf", "kf_converged && good_combinaison && associable && combinaison_well_placed_after_kf");

  std::cout << (float) associable / kfConverged << std::endl;
  std::cout << (float) goodCombinaison / associable << std::endl;
  std::cout << (float) wellPlaced / goodCombinaison << std::endl;
  std::cout << (float) kfConverged / total << std::endl;

}
